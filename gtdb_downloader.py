"""GTDB genome downloader via GTDB metadata + NCBI datasets.

This downloader does NOT fetch the full GTDB representative tarball.
Instead it:

1. Downloads the small GTDB taxonomy and metadata tables for a release.
2. Filters rows by the requested GTDB taxon/taxa.
3. Extracts assembly accession numbers from the metadata table.
4. Splits the accessions into batches.
5. Launches multiple `datasets download genome accession` jobs in parallel.

All network transfers are performed either by Python's stdlib HTTP stack
or by the NCBI `datasets` CLI. No external `wget` / `curl` path is used.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import os
import re
import shutil
import signal
import subprocess
import urllib.error
import urllib.request
import zipfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Dict, Iterable, List, Optional, Set, Tuple

GTDB_BASE = "https://data.gtdb.ecogenomic.org/releases"
DEFAULT_RELEASE = "220.0"

_INTERRUPTED = False


def _sigint(_signum, _frame) -> None:  # pragma: no cover
    global _INTERRUPTED
    _INTERRUPTED = True


def _download_file(url: str, dest: str) -> None:
    """Download a GTDB text file using only Python stdlib HTTP."""
    os.makedirs(os.path.dirname(dest) or ".", exist_ok=True)
    with urllib.request.urlopen(url, timeout=120) as resp, open(dest, "wb") as out:
        shutil.copyfileobj(resp, out, length=4 * 1024 * 1024)


def _http_head_size(url: str, timeout: int = 15) -> Optional[int]:
    req = urllib.request.Request(url, method="HEAD")
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            raw = resp.headers.get("Content-Length")
            return int(raw) if raw and raw.isdigit() else None
    except (urllib.error.URLError, urllib.error.HTTPError, ValueError):
        return None


def _open_text_maybe_gzip(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="ignore")
    return open(path, "r", encoding="utf-8", errors="ignore")


def release_dir_url(release: str) -> str:
    short = release.split(".")[0]
    return f"{GTDB_BASE}/release{short}/{release}/"


def taxonomy_url(release: str, kind: str) -> str:
    return f"{release_dir_url(release)}{kind}_taxonomy_r{release.split('.')[0]}.tsv"


def metadata_url(release: str, kind: str) -> str:
    return f"{release_dir_url(release)}{kind}_metadata_r{release.split('.')[0]}.tsv.gz"


def list_releases() -> List[str]:
    with urllib.request.urlopen(f"{GTDB_BASE}/", timeout=30) as resp:
        html = resp.read().decode("utf-8", errors="ignore")
    rels = sorted({m.group(1) for m in re.finditer(r'href="release(\d+)/"', html)})
    out: List[str] = []
    for rel in rels:
        url = f"{GTDB_BASE}/release{rel}/"
        try:
            with urllib.request.urlopen(url, timeout=15) as resp:
                page = resp.read().decode("utf-8", errors="ignore")
        except (urllib.error.URLError, urllib.error.HTTPError):
            continue
        out.extend(m.group(1) for m in re.finditer(r'href="(\d+\.\d+)/"', page))
    return sorted(set(out))


def _normalize_taxon(text: str) -> str:
    return text.strip().strip("'\"")


def _matches_taxa(taxonomy: str, taxa: List[str], rank_filter: Optional[str]) -> bool:
    parts = [p.strip() for p in taxonomy.split(";") if p.strip()]
    for query in taxa:
        q = _normalize_taxon(query)
        for part in parts:
            if rank_filter:
                prefix = f"{rank_filter}__"
                if not part.startswith(prefix):
                    continue
            if part == q:
                return True
    return False


def download_taxonomy(release: str, output_root: str) -> Tuple[str, str]:
    cache = os.path.join(output_root, "_gtdb_cache", f"r{release.split('.')[0]}", "taxonomy")
    os.makedirs(cache, exist_ok=True)
    paths: List[str] = []
    for kind in ("bac120", "ar53"):
        url = taxonomy_url(release, kind)
        dest = os.path.join(cache, os.path.basename(url))
        if not os.path.exists(dest) or os.path.getsize(dest) == 0:
            print(f"  downloading taxonomy {kind}: {url}", flush=True)
            _download_file(url, dest)
        else:
            print(f"  cached taxonomy: {dest}", flush=True)
        paths.append(dest)
    return tuple(paths)  # type: ignore[return-value]


def download_metadata(release: str, output_root: str) -> Tuple[str, str]:
    cache = os.path.join(output_root, "_gtdb_cache", f"r{release.split('.')[0]}", "metadata")
    os.makedirs(cache, exist_ok=True)
    paths: List[str] = []
    for kind in ("bac120", "ar53"):
        url = metadata_url(release, kind)
        dest = os.path.join(cache, os.path.basename(url))
        if not os.path.exists(dest) or os.path.getsize(dest) < 1024:
            size = _http_head_size(url) or 0
            print(f"  downloading metadata {kind} ({size / 1e6:.1f} MB): {url}", flush=True)
            _download_file(url, dest)
        else:
            print(f"  cached metadata: {dest}", flush=True)
        paths.append(dest)
    return tuple(paths)  # type: ignore[return-value]


def collect_genome_ids(
    taxonomy_paths: Iterable[str],
    taxa: List[str],
    rank: Optional[str] = None,
) -> Set[str]:
    wanted = {_normalize_taxon(t) for t in taxa if t.strip()}
    matches: Set[str] = set()
    for path in taxonomy_paths:
        with _open_text_maybe_gzip(path) as fh:
            for line in fh:
                if not line or line.startswith("accession"):
                    continue
                parts = line.rstrip("\n").split("\t", 1)
                if len(parts) != 2:
                    continue
                genome_id, taxonomy = parts
                if _matches_taxa(taxonomy, list(wanted), rank):
                    matches.add(genome_id.strip())
    return matches


def _resolve_ncbi_accession(header: List[str], row: List[str]) -> str:
    """Extract NCBI assembly accession from GTDB metadata row.

    Priority:
    1. `accession` column stripped from RS_/GB_ prefixes.
    2. `ncbi_genbank_assembly_accession`
    """
    index = {name: i for i, name in enumerate(header)}
    gtdb_acc = row[index["accession"]].strip() if "accession" in index and index["accession"] < len(row) else ""
    if gtdb_acc:
        return re.sub(r"^(RS_|GB_)", "", gtdb_acc)
    if "ncbi_genbank_assembly_accession" in index and index["ncbi_genbank_assembly_accession"] < len(row):
        return row[index["ncbi_genbank_assembly_accession"]].strip()
    return ""


def parse_metadata_to_accession_map(
    metadata_paths: Iterable[str],
    wanted_genome_ids: Set[str],
    only_representative: bool = False,
) -> Dict[str, str]:
    out: Dict[str, str] = {}
    if not wanted_genome_ids:
        return out
    for path in metadata_paths:
        with _open_text_maybe_gzip(path) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            index = {name: i for i, name in enumerate(header)}
            if "accession" not in index:
                continue
            repr_idx = index.get("gtdb_representative", -1)
            for line in fh:
                row = line.rstrip("\n").split("\t")
                if index["accession"] >= len(row):
                    continue
                gtdb_id = row[index["accession"]].strip()
                if gtdb_id not in wanted_genome_ids:
                    continue
                if only_representative and repr_idx >= 0:
                    if repr_idx >= len(row) or row[repr_idx].strip().lower() != "t":
                        continue
                ncbi_acc = _resolve_ncbi_accession(header, row)
                if ncbi_acc:
                    out[gtdb_id] = ncbi_acc
    return out


def _datasets_executable() -> Optional[str]:
    for name in ("datasets", "ncbi-datasets"):
        path = shutil.which(name)
        if path:
            return path
    return None


def _chunked(items: List[str], size: int) -> List[List[str]]:
    return [items[i : i + size] for i in range(0, len(items), size)]


def _run_datasets_batch(
    datasets_bin: str,
    accessions: List[str],
    batch_dir: str,
    batch_name: str,
) -> Optional[str]:
    if not accessions:
        return None
    os.makedirs(batch_dir, exist_ok=True)
    input_txt = os.path.join(batch_dir, f"{batch_name}.txt")
    zip_path = os.path.join(batch_dir, f"{batch_name}.zip")
    log_path = os.path.join(batch_dir, f"{batch_name}.log")
    with open(input_txt, "w") as fh:
        fh.write("\n".join(accessions) + "\n")
    cmd = [
        datasets_bin,
        "download",
        "genome",
        "accession",
        "--inputfile",
        input_txt,
        "--filename",
        zip_path,
        "--no-progressbar",
    ]
    with open(log_path, "w") as logf:
        proc = subprocess.run(cmd, stdout=logf, stderr=subprocess.STDOUT, timeout=1800)
    if proc.returncode != 0:
        return None
    return zip_path if os.path.exists(zip_path) else None


def _extract_datasets_zip(zip_path: str, fasta_dir: str) -> List[str]:
    extracted: List[str] = []
    os.makedirs(fasta_dir, exist_ok=True)
    with zipfile.ZipFile(zip_path, "r") as zf:
        for member in zf.namelist():
            if not member.endswith((".fna", ".fa", ".fasta", ".fna.gz", ".fa.gz", ".fasta.gz")):
                continue
            target = os.path.join(fasta_dir, os.path.basename(member))
            with zf.open(member) as src, open(target, "wb") as dst:
                shutil.copyfileobj(src, dst, length=4 * 1024 * 1024)
            extracted.append(target)
    return extracted


def download_with_datasets(
    accessions: List[str],
    fasta_dir: str,
    threads: int,
    batch_size: int,
) -> Tuple[int, int]:
    datasets_bin = _datasets_executable()
    if not datasets_bin:
        raise RuntimeError("`datasets` CLI not found in the current pixi environment.")
    os.makedirs(fasta_dir, exist_ok=True)
    batch_root = os.path.join(fasta_dir, "_datasets_batches")
    os.makedirs(batch_root, exist_ok=True)

    batches = _chunked(sorted(set(accessions)), max(1, batch_size))
    zip_paths: List[str] = []
    failed_batches = 0

    with ThreadPoolExecutor(max_workers=max(1, threads)) as pool:
        futures = {
            pool.submit(
                _run_datasets_batch,
                datasets_bin,
                batch,
                batch_root,
                f"batch_{i:04d}",
            ): i
            for i, batch in enumerate(batches)
        }
        for future in as_completed(futures):
            if _INTERRUPTED:
                break
            try:
                result = future.result()
            except Exception:
                result = None
            if result:
                zip_paths.append(result)
            else:
                failed_batches += 1

    extracted = 0
    for zip_path in zip_paths:
        extracted += len(_extract_datasets_zip(zip_path, fasta_dir))
    return extracted, failed_batches


def write_accession_list(
    mapping: Dict[str, str],
    output_dir: str,
    release: str,
    query_labels: List[str],
    mode: str,
    source_label: str = "GTDB",
) -> str:
    path = os.path.join(output_dir, "accessions.csv")
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["Release", "Query", "GTDB_Genome_ID", "NCBI_Accession", "Mode", "Source"])
        for gid, acc in sorted(mapping.items()):
            for label in query_labels:
                writer.writerow([release, label, gid, acc, mode, source_label])
    return path


def _sanitize_dir(name: str) -> str:
    name = re.sub(r"[^A-Za-z0-9_.-]+", "_", name.strip())
    return name.strip("_") or "taxon"


def _normalize_accession(text: str) -> str:
    text = text.strip()
    text = re.sub(r"^(RS_|GB_)", "", text)
    return text


def _collect_direct_accessions(raw_accessions: Optional[List[str]], accession_file: Optional[str]) -> List[str]:
    values: List[str] = []
    for item in raw_accessions or []:
        values.extend([x for x in item.split(",") if x.strip()])
    if accession_file:
        with open(accession_file, "r", encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                values.extend([x for x in re.split(r"[\s,]+", line) if x.strip()])
    normalized = [_normalize_accession(x) for x in values if x.strip()]
    return sorted(set(x for x in normalized if x))


def _materialize_import_dir(pool_dir: str, import_to_dir: str, import_mode: str) -> int:
    os.makedirs(import_to_dir, exist_ok=True)
    imported = 0
    for name in os.listdir(pool_dir):
        src = os.path.join(pool_dir, name)
        dst = os.path.join(import_to_dir, name)
        if os.path.isdir(src) or os.path.exists(dst):
            continue
        if import_mode == "copy":
            shutil.copy2(src, dst)
        else:
            try:
                os.symlink(os.path.abspath(src), dst)
            except OSError:
                shutil.copy2(src, dst)
        imported += 1
    return imported


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        prog="gtdb_downloader",
        description=(
            "Resolve GTDB taxon -> GTDB genome IDs -> NCBI accessions from GTDB metadata, "
            "then batch-download genomes with the NCBI `datasets` CLI."
        ),
    )
    parser.add_argument("-R", "--release", default=DEFAULT_RELEASE, help="GTDB release, e.g. 220.0")
    parser.add_argument(
        "-t",
        "--taxon",
        action="append",
        help="GTDB taxon to keep. Repeatable. Comma-separated values inside one --taxon are also accepted.",
    )
    parser.add_argument(
        "-a",
        "--accession",
        action="append",
        help="Direct NCBI assembly accession(s) to download. Repeatable and comma-separated values are supported.",
    )
    parser.add_argument(
        "-A",
        "--accession-file",
        default=None,
        help="Text file containing one accession per line, or comma/space separated accessions.",
    )
    parser.add_argument("-r", "--rank", default=None, choices=["d", "p", "c", "o", "f", "g", "s"])
    parser.add_argument(
        "-m",
        "--scope",
        "--mode",
        choices=["representative", "all", "accessions-only"],
        default="representative",
        help="representative = only GTDB representative genomes; all = all matched genomes; accessions-only = no FASTA download.",
    )
    parser.add_argument("-o", "--out-dir", "--output-root", dest="output_root", default="gtdb_downloads")
    parser.add_argument(
        "-j",
        "--threads",
        type=int,
        default=4,
        help="Number of concurrent `datasets` batch tasks.",
    )
    parser.add_argument(
        "-b",
        "--batch-size",
        type=int,
        default=50,
        help="How many accessions to put into each `datasets download genome accession` task.",
    )
    parser.add_argument(
        "-i",
        "--import-dir",
        "--import-to-dir",
        default=None,
        help="Optional folder to receive the downloaded GTDB fasta files so they can be used directly as input to the dataset-management step.",
    )
    parser.add_argument(
        "--import-mode",
        choices=["symlink", "copy"],
        default="symlink",
        help="How to populate --import-to-dir: symlink is fast and space-saving, copy is fully standalone.",
    )
    parser.add_argument("-L", "--list-releases", action="store_true")
    parser.add_argument("--no-cache", action="store_true")
    args = parser.parse_args(argv)

    if args.list_releases:
        for rel in list_releases():
            print(rel)
        return 0

    if not _datasets_executable() and args.mode != "accessions-only":
        raise SystemExit("`datasets` CLI not found. Please run inside the pixi environment.")

    taxa: List[str] = []
    for item in args.taxon or []:
        taxa.extend([t for t in item.split(",") if t.strip()])
    direct_accessions = _collect_direct_accessions(args.accession, args.accession_file)
    if not taxa and not direct_accessions:
        parser.error("Please provide --taxon or --accession/--accession-file.")

    signal.signal(signal.SIGINT, _sigint)
    started = dt.datetime.now()
    os.makedirs(args.output_root, exist_ok=True)

    if args.no_cache:
        cache_root = os.path.join(args.output_root, "_gtdb_cache")
        if os.path.isdir(cache_root):
            shutil.rmtree(cache_root, ignore_errors=True)

    print("== GTDB genome downloader ==", flush=True)
    print(f"Release: {args.release}", flush=True)
    print(f"Taxa: {taxa}", flush=True)
    print(f"Direct accessions: {len(direct_accessions)}", flush=True)
    print(f"Mode: {args.mode}", flush=True)
    print(f"Threads: {args.threads}", flush=True)
    print(f"Batch size: {args.batch_size}", flush=True)
    print(f"Output root: {os.path.abspath(args.output_root)}", flush=True)
    if args.import_to_dir:
        print(f"Import to dir: {os.path.abspath(args.import_to_dir)} ({args.import_mode})", flush=True)

    genome_ids: Set[str] = set()
    acc_map: Dict[str, str] = {}
    query_labels: List[str] = []
    query_root = "by_taxon"

    if taxa:
        print("\n[1/4] Downloading GTDB taxonomy tables ...", flush=True)
        bac_tx, arc_tx = download_taxonomy(args.release, args.output_root)

        print("\n[2/4] Resolving taxa to GTDB genome IDs ...", flush=True)
        genome_ids = collect_genome_ids([bac_tx, arc_tx], taxa, args.rank)
        print(f"  matched genome IDs: {len(genome_ids)}", flush=True)
        if not genome_ids:
            print("No genomes matched.", flush=True)
            return 1

        print("\n[3/4] Downloading GTDB metadata and extracting assembly accessions ...", flush=True)
        bac_meta, arc_meta = download_metadata(args.release, args.output_root)
        acc_map = parse_metadata_to_accession_map(
            [bac_meta, arc_meta],
            genome_ids,
            only_representative=(args.mode == "representative"),
        )
        print(f"  mapped NCBI accessions: {len(acc_map)}", flush=True)
        if not acc_map:
            print("No NCBI accessions resolved from metadata.", flush=True)
            return 1
        query_labels = taxa
        query_root = "by_taxon"
    else:
        print("\n[1/2] Using direct accession input ...", flush=True)
        acc_map = {acc: acc for acc in direct_accessions}
        query_labels = direct_accessions
        query_root = "by_accession"
        print(f"  normalized direct accessions: {len(acc_map)}", flush=True)

    for label in query_labels:
        query_dir = os.path.join(args.output_root, query_root, _sanitize_dir(label))
        os.makedirs(query_dir, exist_ok=True)
        csv_path = write_accession_list(acc_map, query_dir, args.release, [label], args.mode, source_label="GTDB" if taxa else "Direct_Accession")
        print(f"  wrote accession list -> {csv_path}", flush=True)

    if args.mode == "accessions-only":
        print("\n[4/4] accessions-only mode: FASTA download skipped.", flush=True)
    else:
        print("\n[4/4] Launching batched ncbi_datasets downloads ...", flush=True)
        pool_dir = os.path.join(args.output_root, "fasta_pool")
        extracted, failed_batches = download_with_datasets(
            list(acc_map.values()),
            pool_dir,
            threads=args.threads,
            batch_size=args.batch_size,
        )
        print(f"  extracted FASTA files: {extracted}", flush=True)
        print(f"  failed batches: {failed_batches}", flush=True)
        for label in query_labels:
            query_fasta_dir = os.path.join(args.output_root, query_root, _sanitize_dir(label), "fasta")
            os.makedirs(query_fasta_dir, exist_ok=True)
            for name in os.listdir(pool_dir):
                src = os.path.join(pool_dir, name)
                dst = os.path.join(query_fasta_dir, name)
                if os.path.isdir(src) or os.path.exists(dst):
                    continue
                try:
                    os.symlink(os.path.abspath(src), dst)
                except OSError:
                    shutil.copy2(src, dst)
        if args.import_to_dir:
            imported = _materialize_import_dir(pool_dir, args.import_to_dir, args.import_mode)
            print(f"  imported {imported} fasta files into {args.import_to_dir}", flush=True)

    ended = dt.datetime.now()
    log_path = os.path.join(args.output_root, "gtdb_download.log")
    with open(log_path, "a") as fh:
        fh.write("\n" + "=" * 60 + "\n")
        fh.write(f"GTDB downloader run @ {started.isoformat()}\n")
        fh.write(f"Release: {args.release}\n")
        fh.write(f"Mode: {args.mode}\n")
        fh.write(f"Taxa: {taxa}\n")
        fh.write(f"Direct accessions: {len(direct_accessions)}\n")
        fh.write(f"Threads: {args.threads}\n")
        fh.write(f"Batch size: {args.batch_size}\n")
        fh.write(f"Matched genome IDs: {len(genome_ids)}\n")
        fh.write(f"Mapped NCBI accessions: {len(acc_map)}\n")
        fh.write(f"Import to dir: {args.import_to_dir or ''}\n")
        fh.write(f"Elapsed: {ended - started}\n")
        fh.write("=" * 60 + "\n")
    print(f"\nLog appended to {log_path}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
