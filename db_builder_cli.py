import argparse
import datetime as _dt
import hashlib
import os
import re
import shutil
import subprocess
import unicodedata
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Tuple

import pandas as pd
from tqdm import tqdm


def sanitize_token(text: str, max_len: int, fallback_prefix: str) -> str:
    raw = text if text is not None else ""
    ascii_text = unicodedata.normalize("NFKD", str(raw)).encode("ascii", "ignore").decode("ascii")
    cleaned = re.sub(r"[^A-Za-z0-9]+", "_", ascii_text)
    cleaned = re.sub(r"_+", "_", cleaned).strip("_")
    digest = hashlib.md5(str(raw).encode("utf-8", errors="ignore")).hexdigest()[:8]
    if not cleaned:
        cleaned = f"{fallback_prefix}_{digest}"
    if len(cleaned) > max_len:
        cleaned = f"{cleaned[:max_len-9]}_{digest}"
    return cleaned


def sanitize_alnum_token(text: str, max_len: int, fallback_prefix: str) -> str:
    raw = text if text is not None else ""
    ascii_text = unicodedata.normalize("NFKD", str(raw)).encode("ascii", "ignore").decode("ascii")
    cleaned = re.sub(r"[^A-Za-z0-9]+", "", ascii_text)
    digest = hashlib.md5(str(raw).encode("utf-8", errors="ignore")).hexdigest()[:8].upper()
    fallback_ascii = unicodedata.normalize("NFKD", str(fallback_prefix or "GENOME")).encode("ascii", "ignore").decode("ascii")
    fallback = re.sub(r"[^A-Za-z0-9]+", "", fallback_ascii) or "GENOME"
    if not cleaned:
        cleaned = f"{fallback}{digest}"
    if len(cleaned) > max_len:
        cleaned = cleaned[:max_len]
    return cleaned


def make_db_prefix(db_name: str, max_len: int = 5) -> str:
    token = sanitize_alnum_token(db_name, max_len=max_len, fallback_prefix="BATHY")
    if not token:
        return "Bathy"
    token = token[:max_len]
    return token[:1].upper() + token[1:].lower()


def make_source_letters(source_folder: str, max_len: int = 3) -> str:
    parts = re.findall(r"[A-Za-z]+", str(source_folder or ""))
    initials = "".join(p[0].upper() for p in parts if p)
    if not initials:
        initials = sanitize_alnum_token(source_folder, max_len=max_len, fallback_prefix="SRC").upper()
    if len(initials) < max_len:
        extra = "".join(parts).upper()
        for ch in extra:
            if len(initials) >= max_len:
                break
            if ch not in initials:
                initials += ch
    return initials[:max_len] or "SRC"


def build_renamed_genome_id(db_name: str, source_folder: str, serial: int) -> str:
    db_prefix = make_db_prefix(db_name, max_len=5)
    source_letters = make_source_letters(source_folder, max_len=3)
    return f"{db_prefix}{int(serial):04d}{source_letters}"


def ensure_unique_filename(output_dir: str, base_name: str, ext: str, max_len: int = 120) -> str:
    candidate = f"{base_name}{ext}"
    if not os.path.exists(os.path.join(output_dir, candidate)):
        return candidate
    idx = 2
    while True:
        suffix = f"{idx:02d}"
        trimmed = base_name
        if len(trimmed) + len(suffix) > max_len:
            trimmed = trimmed[: max_len - len(suffix)]
        candidate = f"{trimmed}{suffix}{ext}"
        if not os.path.exists(os.path.join(output_dir, candidate)):
            return candidate
        idx += 1


def has_fasta_extension(filename: str) -> bool:
    name = os.path.basename(filename)
    return name.lower().endswith((".fna", ".fa", ".fasta"))


def is_hidden_file(filename: str) -> bool:
    name = os.path.basename(filename)
    return bool(name) and (name.startswith(".") or name.startswith("._"))


def is_fasta_file(filename: str) -> bool:
    return has_fasta_extension(filename)


def is_probable_text_fasta(path: str, probe_bytes: int = 4096) -> Tuple[bool, str]:
    """Best-effort validation to avoid crashing on AppleDouble/binary junk files."""
    try:
        with open(path, "rb") as fh:
            chunk = fh.read(probe_bytes)
    except OSError as exc:
        return False, f"cannot read file: {exc}"

    if not chunk:
        return False, "empty file"
    if b"\x00" in chunk:
        return False, "binary file detected (NUL byte present)"
    try:
        preview = chunk.decode("utf-8")
    except UnicodeDecodeError:
        return False, "file is not valid UTF-8 text"

    first_nonblank = next((line.strip() for line in preview.splitlines() if line.strip()), "")
    if first_nonblank and not first_nonblank.startswith(">"):
        return False, "first non-empty line is not a FASTA header"
    return True, ""


def classify_input_file(
    path: str,
    exclude_hidden: bool = True,
    strict_fasta_check: bool = True,
) -> Tuple[str, str]:
    name = os.path.basename(path)
    if exclude_hidden and is_hidden_file(name):
        return "hidden_skipped", "hidden file (dotfile or AppleDouble resource fork)"

    has_ext = has_fasta_extension(name)
    if not strict_fasta_check:
        if has_ext:
            return "fasta", ""
        return "non_fasta_ignored", "no standard FASTA suffix and strict FASTA check is disabled"

    looks_like_fasta, reason = is_probable_text_fasta(path)

    if has_ext:
        if strict_fasta_check and not looks_like_fasta:
            return "invalid_fasta_skipped", reason
        return "fasta", ""

    if looks_like_fasta:
        return "fasta_detected_by_content", "no standard FASTA suffix; accepted by content check"

    return "non_fasta_ignored", reason


def get_tool_path(tool_name: str) -> Optional[str]:
    tool_path = shutil.which(tool_name)
    if tool_path:
        return tool_path
    local_tool = os.path.join(os.getcwd(), ".pixi", "bin", tool_name)
    if os.path.exists(local_tool):
        return local_tool
    return None


def resolve_checkm2_db_path() -> Optional[str]:
    """Locate the CheckM2 DIAMOND database. Order:
    1) explicit override passed in (env var CHECKM2DB / DIAMOND_DB, or
       db_path_override argument)
    2) <project>/checkm2_database/CheckM2_database/uniref100.KO.1.dmnd
    3) <project>/checkm2_database/CheckM2_database (directory containing *.dmnd)
    4) <project>/checkm2_database (directory containing *.dmnd)
    Returns the absolute path to the `.dmnd` file (what DIAMOND expects).
    """
    for env_key in ("CHECKM2DB", "DIAMOND_DB"):
        candidate = os.environ.get(env_key)
        if not candidate:
            continue
        if os.path.isfile(candidate):
            return candidate
        if os.path.isdir(candidate):
            dmnd = _find_dmnd_in(candidate)
            if dmnd:
                return dmnd

    project_root = os.path.abspath(os.path.dirname(os.path.abspath(__file__)))
    for candidate in (
        os.path.join(project_root, "checkm2_database", "CheckM2_database"),
        os.path.join(project_root, "checkm2_database"),
    ):
        if os.path.isdir(candidate):
            dmnd = _find_dmnd_in(candidate)
            if dmnd:
                return dmnd
    return None


def _find_dmnd_in(directory: str) -> Optional[str]:
    try:
        for name in os.listdir(directory):
            if name.endswith(".dmnd"):
                return os.path.join(directory, name)
    except OSError:
        return None
    return None


def ensure_checkm2_env() -> Optional[str]:
    db_path = resolve_checkm2_db_path()
    if not db_path:
        return None
    os.environ["CHECKM2DB"] = db_path
    os.environ["DIAMOND_DB"] = db_path
    return db_path


def compute_file_md5(path: str, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.md5()
    with open(path, "rb") as f:
        while True:
            chunk = f.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def compute_fasta_sequence_signature(path: str) -> str:
    digest = hashlib.md5()
    with open(path, "r", encoding="utf-8") as f:
        for raw in f:
            if raw.startswith(">"):
                digest.update(b"\x1e")
                continue
            seq = raw.strip().upper()
            if seq:
                digest.update(seq.encode("ascii", errors="ignore"))
    return digest.hexdigest()


def _numeric_series(df: pd.DataFrame, column: str) -> pd.Series:
    if column not in df.columns:
        return pd.Series(dtype=float)
    return pd.to_numeric(df[column], errors="coerce").dropna()


def _format_numeric_summary(label: str, values: pd.Series, precision: int = 2) -> List[str]:
    if values.empty:
        return [f"{label}: no data"]
    fmt = f"{{:.{precision}f}}"
    return [
        (
            f"{label}: n={len(values)}, min={fmt.format(values.min())}, "
            f"median={fmt.format(values.median())}, mean={fmt.format(values.mean())}, "
            f"max={fmt.format(values.max())}"
        )
    ]


def summarize_metadata_metrics(df: pd.DataFrame) -> Dict[str, List[str]]:
    sections: Dict[str, List[str]] = {"size": [], "quality": [], "rrna16s": []}
    sections["size"].extend(_format_numeric_summary("Genome_Size_bp", _numeric_series(df, "Genome_Size_bp"), precision=0))
    sections["size"].extend(_format_numeric_summary("Contig_Count", _numeric_series(df, "Contig_Count"), precision=0))
    sections["size"].extend(_format_numeric_summary("GC_Content_%", _numeric_series(df, "GC_Content_%"), precision=2))

    completeness = _numeric_series(df, "Completeness")
    contamination = _numeric_series(df, "Contamination")
    sections["quality"].extend(_format_numeric_summary("Completeness", completeness, precision=2))
    if not completeness.empty:
        sections["quality"].append(
            "Completeness bins: "
            f">=90={int((completeness >= 90).sum())}, "
            f"70-90={int(((completeness >= 70) & (completeness < 90)).sum())}, "
            f"50-70={int(((completeness >= 50) & (completeness < 70)).sum())}, "
            f"<50={int((completeness < 50).sum())}"
        )
    sections["quality"].extend(_format_numeric_summary("Contamination", contamination, precision=2))
    if not contamination.empty:
        sections["quality"].append(
            "Contamination bins: "
            f"<5={int((contamination < 5).sum())}, "
            f"5-10={int(((contamination >= 5) & (contamination <= 10)).sum())}, "
            f">10={int((contamination > 10).sum())}"
        )

    rrna_16s = _numeric_series(df, "rRNA_16S_count")
    if rrna_16s.empty:
        sections["rrna16s"].append("16S status: no data")
    else:
        sections["rrna16s"].append(
            "16S status: "
            f"genomes_with_16S={int((rrna_16s > 0).sum())}, "
            f"single_16S={int((rrna_16s == 1).sum())}, "
            f"multi_16S={int((rrna_16s > 1).sum())}, "
            f"without_16S={int((rrna_16s <= 0).sum())}, "
            f"max_copy={int(rrna_16s.max())}"
        )
    return sections


class _UnionFind:
    def __init__(self) -> None:
        self.parent: Dict[str, str] = {}

    def add(self, item: str) -> None:
        if item not in self.parent:
            self.parent[item] = item

    def find(self, item: str) -> str:
        parent = self.parent.setdefault(item, item)
        if parent != item:
            self.parent[item] = self.find(parent)
        return self.parent[item]

    def union(self, a: str, b: str) -> None:
        ra = self.find(a)
        rb = self.find(b)
        if ra != rb:
            self.parent[rb] = ra


def build_exact_duplicate_report(df_metadata: pd.DataFrame, base_out_dir: str, output_csv: str) -> Dict[str, object]:
    signature_to_members: Dict[str, List[str]] = {}
    exact_md5: List[str] = []
    for _, row in df_metadata.iterrows():
        standardized_name = str(row.get("Standardized_Name", "") or "").strip()
        if not standardized_name:
            exact_md5.append("")
            continue
        genome_path = os.path.join(base_out_dir, "genomes", standardized_name)
        if not os.path.exists(genome_path):
            exact_md5.append("")
            continue
        sig = compute_fasta_sequence_signature(genome_path)
        exact_md5.append(sig)
        signature_to_members.setdefault(sig, []).append(standardized_name)
    df_metadata["Exact_Sequence_MD5"] = exact_md5
    df_metadata["Exact_Duplicate_Group"] = ""

    rows: List[Dict[str, object]] = []
    duplicate_group_count = 0
    duplicate_genome_count = 0
    for group_idx, (sig, members) in enumerate(sorted(signature_to_members.items()), start=1):
        if len(members) <= 1:
            continue
        duplicate_group_count += 1
        duplicate_genome_count += len(members)
        group_id = f"ExactDup{group_idx:04d}"
        for standardized_name in members:
            df_metadata.loc[df_metadata["Standardized_Name"] == standardized_name, "Exact_Duplicate_Group"] = group_id
            row = df_metadata.loc[df_metadata["Standardized_Name"] == standardized_name].iloc[0]
            rows.append(
                {
                    "Exact_Duplicate_Group": group_id,
                    "Exact_Sequence_MD5": sig,
                    "Renamed_Genome_ID": row.get("Renamed_Genome_ID", ""),
                    "Standardized_Name": standardized_name,
                    "Original_Name": row.get("Original_Name", ""),
                    "Original_Path": row.get("Original_Path", ""),
                    "Genome_Size_bp": row.get("Genome_Size_bp", ""),
                    "Contig_Count": row.get("Contig_Count", ""),
                    "GC_Content_%": row.get("GC_Content_%", ""),
                }
            )
    pd.DataFrame(rows).to_csv(output_csv, index=False)
    return {
        "output_csv": output_csv,
        "duplicate_group_count": duplicate_group_count,
        "duplicate_genome_count": duplicate_genome_count,
        "unique_signature_count": len(signature_to_members),
    }


def _parse_skani_sparse_value_file(path: str) -> Dict[Tuple[str, str], float]:
    values: Dict[Tuple[str, str], float] = {}
    if not os.path.exists(path):
        return values
    with open(path, "r", encoding="utf-8") as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = re.split(r"\s+", line)
            if len(parts) < 3:
                continue
            a, b = parts[0], parts[1]
            try:
                val = float(parts[-1])
            except ValueError:
                continue
            values[(a, b)] = val
    return values


def run_skani_redundancy_analysis(
    df_metadata: pd.DataFrame,
    base_out_dir: str,
    output_prefix: str,
    threads: int,
    ani_threshold: float,
    af_threshold: float,
) -> Dict[str, object]:
    skani_path = get_tool_path("skani")
    if not skani_path:
        return {"status": "no_tool", "message": "skani executable not found.", "cluster_count": 0}

    genome_paths = [
        os.path.join(base_out_dir, "genomes", str(name))
        for name in df_metadata.get("Standardized_Name", pd.Series(dtype=str)).tolist()
        if str(name).strip()
    ]
    genome_paths = [p for p in genome_paths if os.path.exists(p)]
    if len(genome_paths) < 2:
        return {"status": "not_enough_genomes", "message": "Need at least two genomes for ANI clustering.", "cluster_count": 0}

    output_dir = os.path.dirname(output_prefix)
    os.makedirs(output_dir, exist_ok=True)
    ani_path = output_prefix + ".tsv"
    af_path = ani_path + ".af"
    log_path = output_prefix + ".log"
    list_path = output_prefix + "_genomes.txt"
    with open(list_path, "w", encoding="utf-8") as list_f:
        list_f.write("\n".join(genome_paths) + "\n")
    cmd = [skani_path, "triangle", "-l", list_path, "-E", "-s", "93", "-t", str(int(threads)), "-o", ani_path]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    with open(log_path, "w", encoding="utf-8") as logf:
        logf.write("=== CMD ===\n")
        logf.write(" ".join(cmd) + "\n\n")
        logf.write("=== STDOUT ===\n")
        logf.write(proc.stdout or "")
        logf.write("\n=== STDERR ===\n")
        logf.write(proc.stderr or "")
        logf.write(f"\n=== EXIT CODE ===\n{proc.returncode}\n")
    if proc.returncode != 0:
        return {"status": "failed", "message": f"skani triangle failed; see {log_path}", "cluster_count": 0, "log_path": log_path}
    if not os.path.exists(ani_path):
        return {"status": "no_output", "message": f"skani output missing: {ani_path}", "cluster_count": 0, "log_path": log_path}

    ani_values = _parse_skani_sparse_value_file(ani_path)
    af_values = _parse_skani_sparse_value_file(af_path) if os.path.exists(af_path) else {}
    uf = _UnionFind()
    for standardized_name in df_metadata.get("Standardized_Name", pd.Series(dtype=str)).tolist():
        token = str(standardized_name).strip()
        if token:
            uf.add(token)

    edge_rows: List[Dict[str, object]] = []
    qualifying_edges = 0
    for (a_path, b_path), ani in ani_values.items():
        a_name = os.path.basename(a_path)
        b_name = os.path.basename(b_path)
        af = af_values.get((a_path, b_path), af_values.get((b_path, a_path), float("nan")))
        qualifies = ani >= ani_threshold and (pd.notna(af) and af >= af_threshold)
        edge_rows.append(
            {
                "Genome_A": a_name,
                "Genome_B": b_name,
                "ANI": ani,
                "AF": af if pd.notna(af) else "",
                "ANI_Threshold": ani_threshold,
                "AF_Threshold": af_threshold,
                "Qualified_For_Cluster": bool(qualifies),
            }
        )
        if qualifies:
            qualifying_edges += 1
            uf.union(a_name, b_name)

    edge_csv = output_prefix + "_edges.csv"
    pd.DataFrame(edge_rows).to_csv(edge_csv, index=False)

    cluster_members: Dict[str, List[str]] = {}
    for standardized_name in df_metadata.get("Standardized_Name", pd.Series(dtype=str)).tolist():
        token = str(standardized_name).strip()
        if not token:
            continue
        cluster_members.setdefault(uf.find(token), []).append(token)

    cluster_rows: List[Dict[str, object]] = []
    cluster_lookup: Dict[str, str] = {}
    for cluster_idx, members in enumerate(sorted(cluster_members.values(), key=lambda m: (-len(m), m[0])), start=1):
        cluster_id = f"ANI95AF60_{cluster_idx:04d}"
        for member in sorted(members):
            cluster_lookup[member] = cluster_id
            row = df_metadata.loc[df_metadata["Standardized_Name"] == member].iloc[0]
            cluster_rows.append(
                {
                    "ANI95_AF60_Cluster": cluster_id,
                    "Cluster_Size": len(members),
                    "Renamed_Genome_ID": row.get("Renamed_Genome_ID", ""),
                    "Standardized_Name": member,
                    "Original_Name": row.get("Original_Name", ""),
                    "Original_Path": row.get("Original_Path", ""),
                    "Completeness": row.get("Completeness", ""),
                    "Contamination": row.get("Contamination", ""),
                    "GTDB_Score": row.get("GTDB_Score", ""),
                }
            )
    cluster_csv = output_prefix + "_clusters.csv"
    pd.DataFrame(cluster_rows).to_csv(cluster_csv, index=False)
    df_metadata["ANI95_AF60_Cluster"] = df_metadata["Standardized_Name"].map(cluster_lookup).fillna("")
    multi_member_clusters = sum(1 for members in cluster_members.values() if len(members) > 1)
    clustered_genomes = sum(len(members) for members in cluster_members.values() if len(members) > 1)
    return {
        "status": "ok",
        "message": "",
        "log_path": log_path,
        "ani_path": ani_path,
        "af_path": af_path if os.path.exists(af_path) else "",
        "list_path": list_path,
        "edge_csv": edge_csv,
        "cluster_csv": cluster_csv,
        "cluster_count": len(cluster_members),
        "multi_member_clusters": multi_member_clusters,
        "clustered_genomes": clustered_genomes,
        "qualifying_edges": qualifying_edges,
    }


def standardize_fasta_stream(input_path: str, output_path: str, contig_prefix: str) -> Tuple[int, int, float]:
    contig_idx = 1
    contig_count = 0
    total_len = 0
    gc_count = 0
    with open(input_path, "r", encoding="utf-8") as fin, open(output_path, "w", encoding="utf-8") as fout:
        for line in fin:
            if line.startswith(">"):
                fout.write(f">{contig_prefix}_{contig_idx:04d}\n")
                contig_idx += 1
                contig_count += 1
                continue
            seq = line.strip()
            if not seq:
                continue
            total_len += len(seq)
            s_upper = seq.upper()
            gc_count += s_upper.count("G") + s_upper.count("C")
            fout.write(line)
    if contig_count == 0:
        raise ValueError("no FASTA header found")
    gc_content = (gc_count / total_len * 100.0) if total_len > 0 else 0.0
    return total_len, contig_count, round(gc_content, 2)


def parse_barrnap_gff(gff_path: str) -> Dict[str, int]:
    counts = {"5S": 0, "16S": 0, "23S": 0, "18S": 0, "28S": 0}
    if not os.path.exists(gff_path):
        return {**counts, "total": 0}
    try:
        with open(gff_path, "r") as f:
            for line in f:
                if not line or line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 9:
                    continue
                attrs = parts[8]
                if "5S" in attrs:
                    counts["5S"] += 1
                if "16S" in attrs:
                    counts["16S"] += 1
                if "23S" in attrs:
                    counts["23S"] += 1
                if "18S" in attrs:
                    counts["18S"] += 1
                if "28S" in attrs:
                    counts["28S"] += 1
    except Exception:
        pass
    total = sum(counts.values())
    return {**counts, "total": total}


def _load_fasta_sequences(path: str) -> Dict[str, str]:
    seqs: Dict[str, str] = {}
    cur_id: Optional[str] = None
    cur_chunks: List[str] = []
    with open(path, "r") as f:
        for raw in f:
            line = raw.rstrip("\n").rstrip("\r")
            if not line:
                continue
            if line.startswith(">"):
                if cur_id is not None:
                    seqs[cur_id] = "".join(cur_chunks)
                cur_id = line[1:].split()[0]
                cur_chunks = []
                continue
            cur_chunks.append(line.strip())
    if cur_id is not None:
        seqs[cur_id] = "".join(cur_chunks)
    return seqs


def extract_rrna_sequences(
    fasta_path: str,
    gff_path: str,
    output_fasta: str,
    rrna_type: str = "16S",
    genome_name: str = "",
    source_folder: str = "",
    original_name: str = "",
    original_subfolder: str = "",
) -> int:
    """Slice rRNA sequences out of `fasta_path` using coordinates from `gff_path`.

    Writes records to `output_fasta` and returns the number of records written.
    Works for any rRNA type label that appears in the GFF attributes column.
    """
    if not os.path.exists(gff_path) or not os.path.exists(fasta_path):
        return 0
    seqs = _load_fasta_sequences(fasta_path)
    if not seqs:
        return 0
    target = rrna_type.upper()
    entries: List[str] = []
    os.makedirs(os.path.dirname(output_fasta) or ".", exist_ok=True)
    with open(gff_path, "r") as gff_in:
        for line in gff_in:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            attrs = parts[8]
            if target not in attrs.upper():
                continue
            contig = parts[0]
            if contig not in seqs:
                # contig id from GFF doesn't match fasta; skip.
                continue
            try:
                start = int(parts[3]) - 1  # GFF is 1-based inclusive.
                end = int(parts[4])  # slice end (exclusive).
                strand = parts[6]
            except ValueError:
                continue
            if start < 0:
                start = 0
            if end > len(seqs[contig]):
                end = len(seqs[contig])
            if end <= start:
                continue
            sub = seqs[contig][start:end]
            if strand == "-":
                comp = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
                sub = "".join(comp.get(b.upper(), "N") for b in reversed(sub))
            attrs_dict: Dict[str, str] = {}
            for kv in attrs.split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    attrs_dict[k.strip()] = v.strip()
            entries.append(sub)

    total_entries = len(entries)
    if total_entries == 0:
        return 0

    with open(output_fasta, "w", encoding="utf-8") as out:
        for idx, sub in enumerate(entries, start=1):
            header_id, header_desc = format_rrna_header(
                genome_name=genome_name,
                rrna_type=rrna_type,
                feature_index=idx,
                total_features=total_entries,
                feature_name="",
                contig="",
                start=0,
                end=0,
                strand="",
                product="",
                source_folder=source_folder,
                original_name=original_name,
                original_subfolder=original_subfolder,
            )
            if header_desc:
                out.write(f">{header_id} {header_desc}\n")
            else:
                out.write(f">{header_id}\n")
            for i in range(0, len(sub), 80):
                out.write(sub[i : i + 80] + "\n")
    return total_entries


def run_barrnap(genome_fna: str, output_gff: str, output_rrna_fasta: str, kingdom: str, threads: int) -> Optional[str]:
    barrnap_path = get_tool_path("barrnap")
    if not barrnap_path:
        return None
    os.makedirs(os.path.dirname(output_gff), exist_ok=True)
    log_file = output_gff + ".log"
    cmd = [
        barrnap_path,
        "--kingdom",
        kingdom,
        "--threads",
        str(int(threads)),
        "--outseq",
        output_rrna_fasta,
        genome_fna,
    ]
    try:
        with open(output_gff, "w") as gf, open(log_file, "w") as lf:
            subprocess.run(cmd, check=True, stdout=gf, stderr=lf)
    except subprocess.CalledProcessError:
        return None
    return output_gff


# Project-local CheckM2 1.1.0 + Python 3.12 compatibility shim. The sitecustomize
# inside this directory is auto-imported when its parent is on PYTHONPATH. It
# patches multiprocessing.Pool to use fork and adds unmangled aliases for the
# private Predictor methods so spawn workers can reconstruct them.
CHECKM2_COMPAT_DIR = os.path.join(
    os.path.abspath(os.path.dirname(os.path.abspath(__file__))),
    "checkm2_compat",
)


def run_checkm2(
    genomes_dir: str,
    output_dir: str,
    threads: int,
    status_text: Optional[object] = None,  # reserved for future Streamlit integration
) -> Dict[str, object]:
    """Run CheckM2 `predict` and parse the resulting quality_report.tsv.

    Returns a dict with keys:
      - "status": "ok" | "no_tool" | "no_database" | "failed" | "no_report"
      - "scores":  { <standardized_name.fna>: {completeness, contamination, score} }
      - "db_path": resolved database path or ""
      - "message": human readable message
    """
    checkm2_path = get_tool_path("checkm2")
    if not checkm2_path:
        return {"status": "no_tool", "scores": {}, "db_path": "", "message": "checkm2 executable not found."}

    db_path = ensure_checkm2_env()
    if not db_path:
        msg = (
            "CheckM2 DIAMOND database not found. Expected at "
            "<project>/checkm2_database/CheckM2_database or via $CHECKM2DB / $DIAMOND_DB."
        )
        return {"status": "no_database", "scores": {}, "db_path": "", "message": msg}

    # Surface progress in CLI; ignored when status_text is None (reserved for Streamlit).
    if status_text is not None and hasattr(status_text, "write"):
        try:
            status_text.write("Running CheckM2 predict...")
        except Exception:
            pass

    os.makedirs(output_dir, exist_ok=True)
    log_file = os.path.join(output_dir, "checkm2.log")
    cmd = [
        checkm2_path,
        "predict",
        "--threads",
        str(int(threads)),
        "--input",
        genomes_dir,
        "--output-directory",
        output_dir,
        "-x",
        "fna",
        "--force",
    ]
    # Build a child-process environment that:
    #   * points CheckM2 at the resolved DIAMOND db file
    #   * injects the project-local checkm2_compat/sitecustomize via PYTHONPATH
    child_env = os.environ.copy()
    child_env["CHECKM2DB"] = db_path
    child_env["DIAMOND_DB"] = db_path
    if os.path.isdir(CHECKM2_COMPAT_DIR):
        existing_pp = child_env.get("PYTHONPATH", "")
        child_env["PYTHONPATH"] = (
            CHECKM2_COMPAT_DIR if not existing_pp else CHECKM2_COMPAT_DIR + os.pathsep + existing_pp
        )
    proc = subprocess.run(cmd, capture_output=True, text=True, env=child_env)
    with open(log_file, "w") as f:
        f.write("=== CMD ===\n")
        f.write(" ".join(cmd) + "\n\n")
        f.write("=== STDOUT ===\n")
        f.write(proc.stdout or "")
        f.write("\n=== STDERR ===\n")
        f.write(proc.stderr or "")
        f.write(f"\n=== EXIT CODE ===\n{proc.returncode}\n")
        f.write(f"=== CHECKM2DB (child) ===\n{child_env.get('CHECKM2DB','')}\n")
        f.write(f"=== PYTHONPATH (child) ===\n{child_env.get('PYTHONPATH','')}\n")

    if proc.returncode != 0:
        tail = (proc.stderr or proc.stdout or "").strip().splitlines()
        tail_text = "\n".join(tail[-15:]) if tail else "(no output)"
        return {
            "status": "failed",
            "scores": {},
            "db_path": db_path,
            "message": f"checkm2 predict exited with code {proc.returncode}. Tail:\n{tail_text}",
        }

    report_file = os.path.join(output_dir, "quality_report.tsv")
    if not os.path.exists(report_file):
        return {
            "status": "no_report",
            "scores": {},
            "db_path": db_path,
            "message": f"checkm2 finished but {report_file} is missing.",
        }

    scores, parse_error = parse_checkm2_quality_report(report_file)
    if parse_error:
        return {
            "status": "failed",
            "scores": {},
            "db_path": db_path,
            "message": parse_error,
        }
    return {"status": "ok", "scores": scores, "db_path": db_path, "message": ""}


def parse_checkm2_quality_report(report_file: str) -> Tuple[Dict[str, Dict[str, float]], str]:
    scores: Dict[str, Dict[str, float]] = {}
    try:
        df_q = pd.read_csv(report_file, sep="\t")
        for _, row in df_q.iterrows():
            name = str(row.get("Name", ""))
            if not name.endswith(".fna"):
                name = name + ".fna"
            comp = float(row.get("Completeness", 0))
            contam = float(row.get("Contamination", 100))
            score = comp - 5 * contam
            scores[name] = {"score": score, "completeness": comp, "contamination": contam}
    except Exception as e:
        return {}, f"Failed to parse quality_report.tsv: {e}"
    return scores, ""


def fasta_record_count(path: str) -> int:
    if not os.path.exists(path):
        return 0
    count = 0
    try:
        with open(path, "r", encoding="utf-8") as fh:
            for line in fh:
                if line.startswith(">"):
                    count += 1
    except OSError:
        return 0
    return count


def format_rrna_header(
    genome_name: str,
    rrna_type: str,
    feature_index: int,
    total_features: int,
    feature_name: str,
    contig: str,
    start: int,
    end: int,
    strand: str,
    product: str,
    source_folder: str = "",
    original_name: str = "",
    original_subfolder: str = "",
) -> Tuple[str, str]:
    genome_base = os.path.splitext(os.path.basename(genome_name))[0] if genome_name else "GENOME"
    genome_token = genome_base.strip() or "GENOME"
    del rrna_type, feature_name, contig, start, end, strand, product, source_folder, original_name, original_subfolder
    if int(total_features) <= 1:
        header_id = genome_token
    else:
        header_id = f"{genome_token}_{int(feature_index)}"
    return header_id, ""


def rrna_fasta_has_genome_metadata(path: str, genome_name: str, original_name: str = "") -> bool:
    if not os.path.exists(path):
        return False
    del original_name
    candidate_names = [n for n in (genome_name,) if n]
    candidate_tokens = {
        (os.path.splitext(os.path.basename(name))[0].strip() or "GENOME")
        for name in candidate_names
        if os.path.splitext(os.path.basename(name))[0].strip()
    }
    try:
        with open(path, "r", encoding="utf-8") as fh:
            for line in fh:
                if line.startswith(">"):
                    header = line[1:].strip().split()[0]
                    return bool(header) and (
                        any(header == token or re.fullmatch(rf"{re.escape(token)}_\d+", header) for token in candidate_tokens) or
                        "|16S|" in header or
                        "genome=" in line
                    )
    except OSError:
        return False
    return False


def rewrite_rrna_fasta_headers(
    fasta_path: str,
    record: Dict[str, object],
    rrna_type: str = "16S",
) -> int:
    if not os.path.exists(fasta_path):
        return 0
    tmp_path = fasta_path + ".tmp"
    seq_idx = 0
    genome_name = str(record.get("Standardized_Name", "") or "")
    source_folder = str(record.get("Source_Folder", "") or "")
    original_name = str(record.get("Original_Name", "") or "")
    original_subfolder = str(record.get("Original_Subfolder", "") or "")
    total_records = fasta_record_count(fasta_path)
    try:
        with open(fasta_path, "r", encoding="utf-8") as src, open(tmp_path, "w", encoding="utf-8") as out:
            for line in src:
                if line.startswith(">"):
                    seq_idx += 1
                    raw_header = line[1:].strip()
                    feature_name = raw_header.split()[0] if raw_header else f"{rrna_type}_{seq_idx}"
                    header_id, header_desc = format_rrna_header(
                        genome_name=genome_name,
                        rrna_type=rrna_type,
                        feature_index=seq_idx,
                        total_features=total_records,
                        feature_name=feature_name,
                        contig=str(record.get("Contig_Header_Prefix", "") or genome_name or "CONTIG"),
                        start=seq_idx,
                        end=seq_idx,
                        strand=".",
                        product=raw_header or f"{rrna_type} rRNA",
                        source_folder=source_folder,
                        original_name=original_name,
                        original_subfolder=original_subfolder,
                    )
                    if header_desc:
                        out.write(f">{header_id} {header_desc}\n")
                    else:
                        out.write(f">{header_id}\n")
                else:
                    out.write(line)
        os.replace(tmp_path, fasta_path)
    except OSError:
        try:
            if os.path.exists(tmp_path):
                os.remove(tmp_path)
        except OSError:
            pass
        return 0
    return seq_idx


def append_fasta_file(src_path: str, out_handle) -> int:
    written = 0
    with open(src_path, "r", encoding="utf-8") as src:
        for line in src:
            out_handle.write(line)
            if line.startswith(">"):
                written += 1
    return written


def resolve_record_relpath(record: Dict[str, object], field: str, base_out_dir: str) -> str:
    value = str(record.get(field, "") or "").strip()
    if not value:
        return ""
    path = value if os.path.isabs(value) else os.path.join(base_out_dir, value)
    return path if os.path.exists(path) else ""


def ensure_record_16s_fasta(record: Dict[str, object], base_out_dir: str) -> Tuple[str, int, str]:
    genome_name = str(record.get("Standardized_Name", "") or "").strip()
    contig_prefix = str(record.get("Contig_Header_Prefix", "") or os.path.splitext(genome_name)[0]).strip()
    if not genome_name:
        return "", 0, "missing_standardized_name"

    genome_path = os.path.join(base_out_dir, "genomes", genome_name)
    gff_path = resolve_record_relpath(record, "barrnap_gff", base_out_dir)
    existing_16s = resolve_record_relpath(record, "barrnap_16S_fasta", base_out_dir)
    if existing_16s:
        count = fasta_record_count(existing_16s)
        if count > 0:
            if rrna_fasta_has_genome_metadata(
                existing_16s,
                genome_name,
                original_name=str(record.get("Original_Name", "") or ""),
            ):
                return existing_16s, count, "existing_16s_reused"
            if gff_path and os.path.exists(genome_path):
                # Legacy 16S FASTA lacks genome metadata; regenerate from current GFF/genome
                # so the recovered header retains genome + coordinate information.
                pass
            else:
                rewritten = rewrite_rrna_fasta_headers(existing_16s, record, rrna_type="16S")
                if rewritten > 0:
                    record["rRNA_16S_count"] = rewritten
                    return existing_16s, rewritten, "existing_16s_rewritten"
                return existing_16s, count, "existing_16s_reused"

    if not gff_path or not os.path.exists(genome_path):
        return "", 0, "missing_genome_or_gff"

    out_path = os.path.join(base_out_dir, "barrnap_results", f"{contig_prefix}16S.fasta")
    count = extract_rrna_sequences(
        genome_path,
        gff_path,
        out_path,
        rrna_type="16S",
        genome_name=genome_name,
        source_folder=str(record.get("Source_Folder", "") or ""),
        original_name=str(record.get("Original_Name", "") or ""),
        original_subfolder=str(record.get("Original_Subfolder", "") or ""),
    )
    if count > 0:
        record["barrnap_16S_fasta"] = os.path.relpath(out_path, base_out_dir)
        record["rRNA_16S_count"] = count
        return out_path, count, "generated_from_gff"
    return "", 0, "no_16s_detected"


def collect_database_16s_records(
    records: List[Dict[str, object]],
    base_out_dir: str,
    output_fasta: str,
) -> Dict[str, object]:
    os.makedirs(os.path.dirname(output_fasta) or ".", exist_ok=True)
    genomes_with_16s = 0
    seq_count = 0
    reused_files = 0
    generated_files = 0
    missing_files = 0
    with open(output_fasta, "w", encoding="utf-8") as out:
        for record in records:
            src_path, count, status = ensure_record_16s_fasta(record, base_out_dir)
            if src_path and count > 0:
                seq_count += append_fasta_file(src_path, out)
                genomes_with_16s += 1
                if status == "generated_from_gff":
                    generated_files += 1
                else:
                    reused_files += 1
            else:
                missing_files += 1
    return {
        "output_path": output_fasta,
        "sequence_count": seq_count,
        "genomes_with_16s": genomes_with_16s,
        "reused_files": reused_files,
        "generated_files": generated_files,
        "missing_files": missing_files,
    }


@dataclass
class InputGenome:
    source_folder: str
    source_dir: str
    root: str
    filename: str

    @property
    def path(self) -> str:
        return os.path.join(self.root, self.filename)

    @property
    def subfolder(self) -> str:
        try:
            return os.path.relpath(self.root, self.source_dir)
        except Exception:
            return ""


def iter_input_genomes(source_dirs: List[str]) -> Iterable[InputGenome]:
    for s_dir in source_dirs:
        if not os.path.exists(s_dir):
            continue
        source_label = os.path.basename(os.path.normpath(s_dir))
        for root, _, files in os.walk(s_dir):
            for f in files:
                if not is_fasta_file(f):
                    continue
                yield InputGenome(source_folder=source_label, source_dir=s_dir, root=root, filename=f)


def scan_input_sources(
    source_dirs: List[str],
    exclude_hidden: bool = True,
    strict_fasta_check: bool = True,
) -> Tuple[List[InputGenome], Dict[str, int], int, Dict[str, int], List[Dict[str, str]]]:
    input_genomes: List[InputGenome] = []
    counts: Dict[str, int] = {}
    skipped_records: List[Dict[str, str]] = []
    stats = {
        "hidden_skipped": 0,
        "invalid_fasta_skipped": 0,
        "non_fasta_ignored": 0,
        "fasta_detected_by_content": 0,
    }
    total = 0

    for s_dir in source_dirs:
        if not os.path.exists(s_dir):
            continue
        source_label = os.path.basename(os.path.normpath(s_dir))
        source_count = 0
        for root, _, files in os.walk(s_dir):
            for f in files:
                path = os.path.join(root, f)
                category, reason = classify_input_file(
                    path,
                    exclude_hidden=exclude_hidden,
                    strict_fasta_check=strict_fasta_check,
                )
                if category in {"fasta", "fasta_detected_by_content"}:
                    input_genomes.append(
                        InputGenome(source_folder=source_label, source_dir=s_dir, root=root, filename=f)
                    )
                    source_count += 1
                    total += 1
                    if category == "fasta_detected_by_content":
                        stats["fasta_detected_by_content"] += 1
                    continue
                if category == "hidden_skipped":
                    stats["hidden_skipped"] += 1
                    skipped_records.append(
                        {
                            "Original_Path": path,
                            "Original_Name": f,
                            "Source_Folder": source_label,
                            "Original_Subfolder": os.path.relpath(root, s_dir),
                            "File_MD5": "",
                            "Status": "Hidden_File_Skipped",
                            "Matched_Previous_Name": "",
                            "Reason": reason,
                        }
                    )
                elif category == "invalid_fasta_skipped":
                    stats["invalid_fasta_skipped"] += 1
                    skipped_records.append(
                        {
                            "Original_Path": path,
                            "Original_Name": f,
                            "Source_Folder": source_label,
                            "Original_Subfolder": os.path.relpath(root, s_dir),
                            "File_MD5": "",
                            "Status": "Invalid_FASTA_Skipped",
                            "Matched_Previous_Name": "",
                            "Reason": reason,
                        }
                    )
                else:
                    stats["non_fasta_ignored"] += 1
        counts[source_label] = source_count

    return input_genomes, counts, total, stats, skipped_records


def discover_source_dirs(database_root: str) -> List[str]:
    dirs: List[str] = []
    if not os.path.isdir(database_root):
        return dirs
    for name in sorted(os.listdir(database_root)):
        path = os.path.join(database_root, name)
        if os.path.isdir(path):
            dirs.append(path)
    return dirs


def count_fastas_by_source(source_dirs: List[str]) -> Tuple[Dict[str, int], int]:
    counts: Dict[str, int] = {}
    total = 0
    for s_dir in source_dirs:
        source_label = os.path.basename(os.path.normpath(s_dir))
        count = 0
        for _, _, files in os.walk(s_dir):
            for f in files:
                if is_fasta_file(f):
                    count += 1
        counts[source_label] = count
        total += count
    return counts, total


def find_previous_snapshot(output_root: str, db_name: str, current_version: str) -> Tuple[Optional[str], Optional[str]]:
    if not os.path.isdir(output_root):
        return None, None
    candidates: List[Tuple[float, str, str]] = []
    prefix = f"{db_name}_"
    for name in os.listdir(output_root):
        if not name.startswith(prefix):
            continue
        version = name[len(prefix):]
        if version == current_version:
            continue
        base_dir = os.path.join(output_root, name)
        metadata_path = os.path.join(base_dir, f"{db_name}_{version}_metadata.csv")
        if os.path.exists(metadata_path):
            candidates.append((os.path.getmtime(metadata_path), base_dir, metadata_path))
    if not candidates:
        return None, None
    candidates.sort(key=lambda x: x[0], reverse=True)
    _, base_dir, metadata_path = candidates[0]
    return base_dir, metadata_path


def write_summary_log(path: str, lines: List[str]) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write("\n".join(lines))


def main():
    parser = argparse.ArgumentParser(
        prog="db_builder_cli",
        description="Build a versioned genome database from multiple source folders without modifying original files.",
    )
    parser.add_argument("-n", "--name", "--db-name", dest="db_name", required=True)
    parser.add_argument("-v", "--version", "--db-version", dest="db_version", required=True)
    parser.add_argument("-D", "--source-root", "--database-root", dest="database_root", help="Root directory containing one subfolder per source.")
    parser.add_argument("-s", "--source-dir", "--sources", dest="sources", nargs="+", help="One or more source directories.")
    parser.add_argument("-o", "--out-dir", "--output-root", dest="output_root", default="local_databases", help="Root output directory.")
    parser.add_argument(
        "-j",
        "--threads",
        type=int,
        default=8,
        help="Threads used by barrnap and CheckM2 (barrnap scales well; CheckM2 1.1.0 has a Python 3.12 pickling bug at high thread counts, so 4-8 is a safe ceiling).",
    )
    parser.add_argument("-Q", "--checkm2", "--run-checkm2", dest="run_checkm2", action="store_true", default=False)
    parser.add_argument("-B", "--barrnap", "--run-barrnap", dest="run_barrnap", action="store_true", default=False)
    parser.add_argument("-k", "--rrna-kingdom", "--barrnap-kingdom", dest="barrnap_kingdom", default="bac", choices=["bac", "arc", "euk", "mito"])
    parser.add_argument(
        "-c",
        "--checkm2-db",
        "--checkm2-db-path",
        dest="checkm2_db_path",
        default=None,
        help=(
            "Explicit path to the CheckM2 DIAMOND database file (*.dmnd) or to a "
            "directory containing one. Overrides $CHECKM2DB and the project-local "
            "<project>/checkm2_database/ lookup."
        ),
    )
    parser.add_argument(
        "--exclude-hidden",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Whether to skip hidden files such as .DS_Store and AppleDouble ._* files. Default: enabled.",
    )
    parser.add_argument(
        "--strict-fasta-check",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Validate FASTA format from file content during discovery, and also detect FASTA files without standard suffixes. Default: enabled.",
    )
    parser.add_argument(
        "-S",
        "--collect-16s-to",
        dest="collect_16s_to",
        default=None,
        help="Collect all database-wide 16S sequences into one FASTA file. Existing per-genome 16S results are reused when possible.",
    )
    parser.add_argument(
        "--ani-cluster",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Run skani-based ANI redundancy analysis and report 95%% ANI / AF 60%% clusters. Default: disabled.",
    )
    parser.add_argument("--ani-threshold", type=float, default=95.0, help="ANI threshold for clustering. Default: 95.0.")
    parser.add_argument("--af-threshold", type=float, default=60.0, help="AF threshold for clustering. Default: 60.0.")
    args = parser.parse_args()

    start = _dt.datetime.now()
    base_out_dir = os.path.join(args.output_root, f"{args.db_name}_{args.db_version}")
    genomes_out_dir = os.path.join(base_out_dir, "genomes")
    checkm2_out_dir = os.path.join(base_out_dir, "checkm2_results")
    barrnap_out_dir = os.path.join(base_out_dir, "barrnap_results")
    new_genomes_out_dir = os.path.join(base_out_dir, "_new_genomes_for_qc")
    redundancy_dir = os.path.join(base_out_dir, "redundancy_analysis")
    metadata_path = os.path.join(base_out_dir, f"{args.db_name}_{args.db_version}_metadata.csv")
    mapping_path = os.path.join(base_out_dir, f"{args.db_name}_{args.db_version}_genome_id_mapping.csv")
    scan_path = os.path.join(base_out_dir, f"{args.db_name}_{args.db_version}_scan_inventory.csv")
    source_inventory_path = os.path.join(base_out_dir, f"{args.db_name}_{args.db_version}_source_counts.csv")
    version_compare_path = os.path.join(base_out_dir, f"{args.db_name}_{args.db_version}_version_comparison.csv")
    exact_duplicate_path = os.path.join(redundancy_dir, f"{args.db_name}_{args.db_version}_exact_duplicates.csv")
    ani_output_prefix = os.path.join(redundancy_dir, f"{args.db_name}_{args.db_version}_ani95_af60")
    os.makedirs(genomes_out_dir, exist_ok=True)
    os.makedirs(redundancy_dir, exist_ok=True)
    if args.run_barrnap:
        os.makedirs(barrnap_out_dir, exist_ok=True)

    if args.database_root:
        source_dirs = discover_source_dirs(args.database_root)
    else:
        source_dirs = args.sources or []
    source_counts: Dict[str, int] = {}
    total_input_found = 0
    discovery_stats = {
        "hidden_skipped": 0,
        "invalid_fasta_skipped": 0,
        "non_fasta_ignored": 0,
        "fasta_detected_by_content": 0,
    }
    discovery_skips: List[Dict[str, str]] = []
    input_genomes: List[InputGenome] = []
    current_df = pd.DataFrame()
    current_md5_map: Dict[str, Dict[str, object]] = {}
    resumed_current_version = 0

    if source_dirs:
        input_genomes, source_counts, total_input_found, discovery_stats, discovery_skips = scan_input_sources(
            source_dirs,
            exclude_hidden=args.exclude_hidden,
            strict_fasta_check=args.strict_fasta_check,
        )
        print("== Source FASTA Inventory ==")
        for source_name, count in source_counts.items():
            print(f"  {source_name}: {count}")
        print(f"  Total usable FASTA: {total_input_found}")
        print("== Discovery Summary ==")
        print(
            "  Hidden/invalid FASTA skipped: "
            f"{discovery_stats['hidden_skipped'] + discovery_stats['invalid_fasta_skipped']}"
        )
        print(f"  Hidden files skipped: {discovery_stats['hidden_skipped']}")
        print(f"  Invalid FASTA skipped: {discovery_stats['invalid_fasta_skipped']}")
        print(f"  FASTA detected by content: {discovery_stats['fasta_detected_by_content']}")
        print(f"  Non-FASTA files ignored: {discovery_stats['non_fasta_ignored']}")
        print(f"  Exclude hidden files: {args.exclude_hidden}")
        print(f"  Strict FASTA check: {args.strict_fasta_check}")
    elif os.path.exists(metadata_path):
        print("No source directories provided; reusing existing current-version metadata for resume-only operations.")
    else:
        raise SystemExit("Please provide --database-root or --sources, or reuse an existing current version with --collect-16s-to.")

    prev_base_dir, prev_metadata_path = find_previous_snapshot(args.output_root, args.db_name, args.db_version)
    previous_df = pd.DataFrame()
    previous_md5_map: Dict[str, Dict[str, object]] = {}
    carried_forward = 0
    if os.path.exists(metadata_path):
        current_df = pd.read_csv(metadata_path)
        if "File_MD5" in current_df.columns:
            current_md5_map = {
                str(row["File_MD5"]): row.to_dict()
                for _, row in current_df.dropna(subset=["File_MD5"]).iterrows()
            }
        resumed_current_version = len(current_df)
        print(f"Found current version snapshot: {base_out_dir}")
    elif prev_metadata_path and os.path.exists(prev_metadata_path):
        previous_df = pd.read_csv(prev_metadata_path)
        if "File_MD5" in previous_df.columns:
            previous_md5_map = {
                str(row["File_MD5"]): row.to_dict()
                for _, row in previous_df.dropna(subset=["File_MD5"]).iterrows()
            }
        print(f"Found previous version snapshot: {prev_base_dir}")
        prev_genomes_dir = os.path.join(prev_base_dir, "genomes")
        if os.path.isdir(prev_genomes_dir):
            shutil.copytree(prev_genomes_dir, genomes_out_dir, dirs_exist_ok=True)
        prev_barrnap_dir = os.path.join(prev_base_dir, "barrnap_results")
        if os.path.isdir(prev_barrnap_dir):
            shutil.copytree(prev_barrnap_dir, barrnap_out_dir, dirs_exist_ok=True)
        if not previous_df.empty:
            previous_df["Database_Version"] = args.db_version
            previous_df["Version_Status"] = "Carried_Forward"
            carried_forward = len(previous_df)
    else:
        print("No previous version snapshot found. Building from scratch.")

    if not source_dirs and current_df.empty:
        raise SystemExit("No source data or resumable metadata found.")
    if source_dirs and not input_genomes and current_df.empty:
        raise SystemExit("No fasta files found in the provided source folders.")

    records: List[Dict[str, object]]
    if not current_df.empty:
        records = current_df.to_dict("records")
    else:
        records = previous_df.to_dict("records") if not previous_df.empty else []
    for r in records:
        # Mark carried-forward rows so a missing score is never mistaken for a real N/A.
        r.setdefault("Quality_Status", "carried_forward")
        r.setdefault("Quality_Message", "Inherited from previous version snapshot; not re-evaluated.")
    scan_records: List[Dict[str, object]] = list(discovery_skips)
    new_records: List[Dict[str, object]] = []
    os.makedirs(new_genomes_out_dir, exist_ok=True)
    next_genome_serial = len(records) + 1

    total_items = len(input_genomes)
    progress_counts = {"new": 0, "current": 0, "existing": 0, "invalid": 0, "unreadable": 0}
    progress_bar = tqdm(total=total_items, desc="Genome processing", unit="genome", dynamic_ncols=True) if total_items else None

    def mark_progress(kind: str) -> None:
        progress_counts[kind] = progress_counts.get(kind, 0) + 1
        if progress_bar is not None:
            progress_bar.update(1)
            progress_bar.set_postfix(
                new=progress_counts["new"],
                resumed=progress_counts["current"],
                carried=progress_counts["existing"],
                skipped=progress_counts["invalid"] + progress_counts["unreadable"],
                refresh=False,
            )

    for idx, item in enumerate(input_genomes, start=1):
        if args.strict_fasta_check:
            is_valid_text, invalid_reason = is_probable_text_fasta(item.path)
            if not is_valid_text:
                mark_progress("invalid")
                scan_records.append(
                    {
                        "Original_Path": item.path,
                        "Original_Name": item.filename,
                        "Source_Folder": item.source_folder,
                        "Original_Subfolder": item.subfolder,
                        "File_MD5": "",
                        "Status": "Invalid_FASTA_Skipped",
                        "Matched_Previous_Name": "",
                        "Reason": invalid_reason,
                    }
                )
                continue

        raw_md5 = compute_file_md5(item.path)
        if raw_md5 in current_md5_map:
            cur_row = current_md5_map[raw_md5]
            mark_progress("current")
            scan_records.append(
                {
                    "Original_Path": item.path,
                    "Original_Name": item.filename,
                    "Source_Folder": item.source_folder,
                    "Original_Subfolder": item.subfolder,
                    "File_MD5": raw_md5,
                    "Status": "Current_Version_Reused",
                    "Matched_Previous_Name": cur_row.get("Standardized_Name", ""),
                }
            )
            continue
        if raw_md5 in previous_md5_map:
            prev_row = previous_md5_map[raw_md5]
            mark_progress("existing")
            scan_records.append(
                {
                    "Original_Path": item.path,
                    "Original_Name": item.filename,
                    "Source_Folder": item.source_folder,
                    "Original_Subfolder": item.subfolder,
                    "File_MD5": raw_md5,
                    "Status": "Existing_Skipped",
                    "Matched_Previous_Name": prev_row.get("Standardized_Name", ""),
                }
            )
            continue

        source_code = make_source_letters(item.source_folder, max_len=3)
        renamed_genome_id = build_renamed_genome_id(args.db_name, item.source_folder, next_genome_serial)
        next_genome_serial += 1
        new_filename = ensure_unique_filename(genomes_out_dir, renamed_genome_id, ".fna")
        out_path = os.path.join(genomes_out_dir, new_filename)
        contig_prefix = os.path.splitext(new_filename)[0]
        qc_stage_path = os.path.join(new_genomes_out_dir, new_filename)

        try:
            size_bp, contigs, gc = standardize_fasta_stream(item.path, out_path, contig_prefix)
        except (UnicodeDecodeError, ValueError, OSError) as exc:
            mark_progress("unreadable")
            try:
                if os.path.exists(out_path):
                    os.remove(out_path)
            except OSError:
                pass
            scan_records.append(
                {
                    "Original_Path": item.path,
                    "Original_Name": item.filename,
                    "Source_Folder": item.source_folder,
                    "Original_Subfolder": item.subfolder,
                    "File_MD5": raw_md5,
                    "Status": "Unreadable_FASTA_Skipped",
                    "Matched_Previous_Name": "",
                    "Reason": str(exc),
                }
            )
            continue
        shutil.copy2(out_path, qc_stage_path)

        rrna_counts = {"5S": 0, "16S": 0, "23S": 0, "18S": 0, "28S": 0, "total": 0}
        rrna_gff = ""
        rrna_fasta_rel = ""
        rrna_16s_rel = ""
        if args.run_barrnap:
            rrna_gff_path = os.path.join(barrnap_out_dir, f"{contig_prefix}Barrnap.gff")
            rrna_fasta_path = os.path.join(barrnap_out_dir, f"{contig_prefix}Rrna.fasta")
            rrna_16s_path = os.path.join(barrnap_out_dir, f"{contig_prefix}16S.fasta")
            if os.path.exists(rrna_gff_path) and os.path.exists(rrna_fasta_path):
                print(f"[{idx}/{total_items}] REUSE barrnap results: {rrna_gff_path}")
                gff_out = rrna_gff_path
            else:
                gff_out = run_barrnap(out_path, rrna_gff_path, rrna_fasta_path, args.barrnap_kingdom, args.threads)
            if gff_out:
                rrna_counts = parse_barrnap_gff(gff_out)
                rrna_gff = os.path.relpath(gff_out, base_out_dir)
                rrna_fasta_rel = os.path.relpath(rrna_fasta_path, base_out_dir) if os.path.exists(rrna_fasta_path) else ""
                # Slice 16S rRNA sequences directly from the genome fasta using the
                # coordinates barrnap predicted in the GFF. This is independent of
                # barrnap's own --outseq which only emits the full rRNA mix.
                if os.path.exists(rrna_16s_path) and fasta_record_count(rrna_16s_path) > 0:
                    n_16s = fasta_record_count(rrna_16s_path)
                else:
                    n_16s = extract_rrna_sequences(
                        out_path,
                        gff_out,
                        rrna_16s_path,
                        rrna_type="16S",
                        genome_name=contig_prefix,
                        source_folder=item.source_folder,
                        original_name=item.filename,
                        original_subfolder=item.subfolder,
                    )
                if n_16s > 0:
                    rrna_16s_rel = os.path.relpath(rrna_16s_path, base_out_dir)
                    rrna_counts["16S"] = n_16s

        record = {
            "Database_Version": args.db_version,
            "Renamed_Genome_ID": contig_prefix,
            "Standardized_Name": new_filename,
            "Original_Name": item.filename,
            "Original_Path": item.path,
            "Original_Subfolder": item.subfolder,
            "Source_Folder": item.source_folder,
            "Source_Code": source_code,
            "Genome_Serial": next_genome_serial - 1,
            "Sequence_Header_Prefix": contig_prefix,
            "Contig_Header_Prefix": contig_prefix,
            "File_MD5": raw_md5,
            "Genome_Size_bp": size_bp,
            "Contig_Count": contigs,
            "GC_Content_%": gc,
            "rRNA_5S_count": rrna_counts.get("5S", 0),
            "rRNA_16S_count": rrna_counts.get("16S", 0),
            "rRNA_23S_count": rrna_counts.get("23S", 0),
            "rRNA_18S_count": rrna_counts.get("18S", 0),
            "rRNA_28S_count": rrna_counts.get("28S", 0),
            "rRNA_total": rrna_counts.get("total", 0),
            "barrnap_gff": rrna_gff,
            "barrnap_rrna_fasta": rrna_fasta_rel,
            "barrnap_16S_fasta": rrna_16s_rel,
            "Version_Status": "New_In_This_Version",
        }
        new_records.append(record)
        scan_records.append(
            {
                "Original_Path": item.path,
                "Original_Name": item.filename,
                "Source_Folder": item.source_folder,
                "Original_Subfolder": item.subfolder,
                "File_MD5": raw_md5,
                "Status": "New_Processed",
                "Matched_Previous_Name": "",
                "Renamed_Genome_ID": contig_prefix,
                "Standardized_Name": new_filename,
            }
        )
        mark_progress("new")

    if progress_bar is not None:
        progress_bar.close()
        print(
            "Progress summary: "
            f"new={progress_counts['new']}, "
            f"current_resumed={progress_counts['current']}, "
            f"carried_forward_hits={progress_counts['existing']}, "
            f"invalid_or_unreadable={progress_counts['invalid'] + progress_counts['unreadable']}"
        )

    checkm2_result: Dict[str, object] = {"status": "skipped", "scores": {}, "db_path": "", "message": ""}
    if args.run_checkm2:
        if not new_records:
            checkm2_result = {
                "status": "skipped",
                "scores": {},
                "db_path": "",
                "message": "No new genomes to evaluate (incremental build reused prior/current results).",
            }
            print("CheckM2: skipped (no new genomes in this version).")
        else:
            # If user passed --checkm2-db-path explicitly, set env vars up front.
            if args.checkm2_db_path:
                p = os.path.abspath(os.path.expanduser(args.checkm2_db_path))
                if not (os.path.isfile(p) or os.path.isdir(p)):
                    print(
                        f"WARNING: --checkm2-db-path {p} does not exist; "
                        "falling back to auto-discovery."
                    )
                else:
                    os.environ["CHECKM2DB"] = p
                    os.environ["DIAMOND_DB"] = p
            pre_db = ensure_checkm2_env()
            if pre_db:
                print(f"CheckM2 database located: {pre_db}")
            else:
                print(
                    "CheckM2 database NOT located. Searched: "
                    "--checkm2-db-path, <project>/checkm2_database/CheckM2_database, "
                    "$CHECKM2DB, $DIAMOND_DB."
                )
            existing_report = os.path.join(checkm2_out_dir, "quality_report.tsv")
            expected_names = {str(r["Standardized_Name"]) for r in new_records}
            existing_scores: Dict[str, Dict[str, float]] = {}
            if os.path.exists(existing_report):
                existing_scores, parse_error = parse_checkm2_quality_report(existing_report)
                if parse_error:
                    print(f"Existing CheckM2 report ignored: {parse_error}")
                    existing_scores = {}
            if expected_names and expected_names.issubset(existing_scores.keys()):
                print(f"REUSE existing CheckM2 report: {existing_report}")
                checkm2_result = {
                    "status": "resumed",
                    "scores": {k: existing_scores[k] for k in expected_names},
                    "db_path": pre_db or "",
                    "message": "Existing CheckM2 report already covers all genomes staged for this run.",
                }
            else:
                print(f"Running CheckM2 on {len(new_records)} new genomes...")
                checkm2_result = run_checkm2(new_genomes_out_dir, checkm2_out_dir, args.threads)
            print(f"CheckM2 status: {checkm2_result.get('status')}")
            if checkm2_result.get("message"):
                print(f"CheckM2 message: {checkm2_result['message']}")

    checkm2_scores = checkm2_result.get("scores", {}) if isinstance(checkm2_result, dict) else {}
    quality_status = str(checkm2_result.get("status", "skipped")) if isinstance(checkm2_result, dict) else "skipped"
    quality_message = str(checkm2_result.get("message", "")) if isinstance(checkm2_result, dict) else ""
    if not checkm2_scores and args.run_checkm2 and new_records:
        # Make the failure visible in the CSV.
        for r in new_records:
            r.setdefault("Completeness", "N/A")
            r.setdefault("Contamination", "N/A")
            r.setdefault("GTDB_Score", "N/A")

    for r in new_records:
        g_name = str(r["Standardized_Name"])
        c_info = checkm2_scores.get(g_name, {}) if isinstance(checkm2_scores, dict) else {}
        if isinstance(c_info, dict) and c_info:
            r["Completeness"] = c_info.get("completeness", "N/A")
            r["Contamination"] = c_info.get("contamination", "N/A")
            r["GTDB_Score"] = c_info.get("score", "N/A")
            r["Quality_Status"] = "ok"
            r["Quality_Message"] = ""
        else:
            r["Completeness"] = "N/A"
            r["Contamination"] = "N/A"
            r["GTDB_Score"] = "N/A"
            r["Quality_Status"] = quality_status if args.run_checkm2 else "not_run"
            r["Quality_Message"] = quality_message if args.run_checkm2 else "CheckM2 not requested."

    records.extend(new_records)
    df_metadata = pd.DataFrame(records)
    for col in [
        "Renamed_Genome_ID",
        "Source_Code",
        "Genome_Serial",
        "Sequence_Header_Prefix",
        "Original_Name",
        "Original_Path",
        "Original_Subfolder",
        "Source_Folder",
        "Standardized_Name",
        "Contig_Header_Prefix",
        "File_MD5",
        "rRNA_total",
        "Version_Status",
        "Database_Version",
        "Quality_Status",
        "Quality_Message",
        "barrnap_16S_fasta",
        "Genome_Size_bp",
        "Contig_Count",
        "GC_Content_%",
        "Completeness",
        "Contamination",
        "GTDB_Score",
        "Exact_Sequence_MD5",
        "Exact_Duplicate_Group",
        "ANI95_AF60_Cluster",
    ]:
        if col not in df_metadata.columns:
            df_metadata[col] = ""

    scan_df = pd.DataFrame(scan_records)
    scan_df.to_csv(scan_path, index=False)

    source_inventory_df = pd.DataFrame(
        [{"Source_Folder": source, "Fasta_Count": count} for source, count in source_counts.items()]
    )
    source_inventory_df.to_csv(source_inventory_path, index=False)

    version_compare_rows = []
    previous_total = carried_forward
    current_total = len(df_metadata)
    version_compare_rows.append(
        {
            "DB_Name": args.db_name,
            "Current_Version": args.db_version,
            "Previous_Version_Path": prev_base_dir or "",
            "Previous_Total_Genomes": previous_total,
            "New_Genomes_Processed": len(new_records),
            "Existing_Genomes_Skipped": int((scan_df["Status"] == "Existing_Skipped").sum()) if not scan_df.empty else 0,
            "Current_Total_Genomes": current_total,
        }
    )
    pd.DataFrame(version_compare_rows).to_csv(version_compare_path, index=False)

    collect_16s_result: Optional[Dict[str, object]] = None
    if args.collect_16s_to:
        collect_path = os.path.abspath(os.path.expanduser(args.collect_16s_to))
        collect_16s_result = collect_database_16s_records(records, base_out_dir, collect_path)
        df_metadata = pd.DataFrame(records)
        print(
            "Collected database-wide 16S: "
            f"{collect_16s_result['sequence_count']} sequences from "
            f"{collect_16s_result['genomes_with_16s']} genomes -> {collect_path}"
        )

    exact_duplicate_result = build_exact_duplicate_report(df_metadata, base_out_dir, exact_duplicate_path)
    ani_cluster_result: Dict[str, object] = {
        "status": "not_requested",
        "message": "ANI clustering not requested.",
        "cluster_count": 0,
    }
    if args.ani_cluster:
        print(
            f"Running ANI clustering with skani (ANI>={args.ani_threshold}, AF>={args.af_threshold})..."
        )
        ani_cluster_result = run_skani_redundancy_analysis(
            df_metadata,
            base_out_dir,
            ani_output_prefix,
            args.threads,
            args.ani_threshold,
            args.af_threshold,
        )
        print(f"ANI clustering status: {ani_cluster_result.get('status')}")
        if ani_cluster_result.get("message"):
            print(f"ANI clustering message: {ani_cluster_result['message']}")

    df_metadata.to_csv(metadata_path, index=False)

    mapping_df = df_metadata[
        [
            "Renamed_Genome_ID",
            "Standardized_Name",
            "Source_Code",
            "Genome_Serial",
            "Sequence_Header_Prefix",
            "Contig_Header_Prefix",
            "Original_Name",
            "Original_Path",
            "Original_Subfolder",
            "Source_Folder",
            "File_MD5",
            "Genome_Size_bp",
            "Contig_Count",
            "GC_Content_%",
            "Completeness",
            "Contamination",
            "GTDB_Score",
            "rRNA_total",
            "rRNA_16S_count",
            "Exact_Duplicate_Group",
            "ANI95_AF60_Cluster",
            "Version_Status",
            "Database_Version",
        ]
    ].copy()
    mapping_df.to_csv(mapping_path, index=False)

    metric_sections = summarize_metadata_metrics(df_metadata)

    end = _dt.datetime.now()
    summary_lines = []
    summary_lines.append("=" * 60)
    summary_lines.append("GTDB Renew - Database Builder (CLI)")
    summary_lines.append("=" * 60)
    summary_lines.append(f"DB Name: {args.db_name}")
    summary_lines.append(f"DB Version: {args.db_version}")
    summary_lines.append(f"Start Time: {start.strftime('%Y-%m-%d %H:%M:%S')}")
    summary_lines.append(f"End Time: {end.strftime('%Y-%m-%d %H:%M:%S')}")
    summary_lines.append(f"Elapsed: {str(end - start)}")
    summary_lines.append("")
    summary_lines.append("--- Parameters ---")
    summary_lines.append(f"Threads: {args.threads}")
    summary_lines.append(f"Run CheckM2: {args.run_checkm2}")
    summary_lines.append(f"Run barrnap: {args.run_barrnap}")
    summary_lines.append(f"Collect all 16S to: {args.collect_16s_to or ''}")
    summary_lines.append(f"Run ANI clustering: {args.ani_cluster}")
    summary_lines.append(f"ANI threshold: {args.ani_threshold}")
    summary_lines.append(f"AF threshold: {args.af_threshold}")
    summary_lines.append(f"Exclude hidden files: {args.exclude_hidden}")
    summary_lines.append(f"Strict FASTA check: {args.strict_fasta_check}")
    if args.run_barrnap:
        summary_lines.append(f"barrnap kingdom: {args.barrnap_kingdom}")
    summary_lines.append(f"Database root: {args.database_root or ''}")
    if args.run_checkm2:
        summary_lines.append("")
        summary_lines.append("--- CheckM2 ---")
        summary_lines.append(f"Status: {checkm2_result.get('status', 'skipped')}")
        db_path = checkm2_result.get('db_path', '')
        summary_lines.append(f"Database path: {db_path or '(NOT FOUND)'}")
        msg = checkm2_result.get('message', '')
        if msg:
            summary_lines.append("Message:")
            for line in str(msg).splitlines():
                summary_lines.append(f"  {line}")
        else:
            summary_lines.append(f"Genomes evaluated: {len(checkm2_result.get('scores', {}))}")
        summary_lines.append(f"CheckM2 log: {os.path.abspath(os.path.join(checkm2_out_dir, 'checkm2.log'))}")
    summary_lines.append("")
    summary_lines.append("--- Input ---")
    summary_lines.append(f"Source folders: {len(source_dirs)}")
    summary_lines.append(f"Input genomes found: {len(input_genomes)}")
    summary_lines.append(
        "Hidden/invalid FASTA skipped: "
        f"{discovery_stats['hidden_skipped'] + discovery_stats['invalid_fasta_skipped']}"
    )
    summary_lines.append(f"Hidden files skipped: {discovery_stats['hidden_skipped']}")
    summary_lines.append(f"Invalid FASTA skipped: {discovery_stats['invalid_fasta_skipped']}")
    summary_lines.append(f"FASTA detected by content: {discovery_stats['fasta_detected_by_content']}")
    summary_lines.append(f"Non-FASTA files ignored: {discovery_stats['non_fasta_ignored']}")
    for source_name, count in source_counts.items():
        summary_lines.append(f"  - {source_name}: {count}")
    summary_lines.append("")
    summary_lines.append("--- Incremental Mode ---")
    summary_lines.append(f"Previous snapshot: {prev_base_dir or 'None'}")
    summary_lines.append(f"Current version rows resumed: {resumed_current_version}")
    summary_lines.append(f"Carried forward genomes: {carried_forward}")
    summary_lines.append(f"New genomes processed: {len(new_records)}")
    summary_lines.append(f"Existing genomes skipped: {int((scan_df['Status'] == 'Existing_Skipped').sum()) if not scan_df.empty else 0}")
    summary_lines.append(f"Current version reused: {int((scan_df['Status'] == 'Current_Version_Reused').sum()) if not scan_df.empty else 0}")
    if collect_16s_result is not None:
        summary_lines.append("")
        summary_lines.append("--- 16S Collection ---")
        summary_lines.append(f"Output FASTA: {collect_16s_result['output_path']}")
        summary_lines.append(f"Sequences written: {collect_16s_result['sequence_count']}")
        summary_lines.append(f"Genomes with 16S: {collect_16s_result['genomes_with_16s']}")
        summary_lines.append(f"Existing 16S reused: {collect_16s_result['reused_files']}")
        summary_lines.append(f"16S regenerated from GFF: {collect_16s_result['generated_files']}")
        summary_lines.append(f"Genomes without 16S: {collect_16s_result['missing_files']}")
    summary_lines.append("")
    summary_lines.append("--- Sequence Metrics ---")
    for line in metric_sections["size"]:
        summary_lines.append(line)
    summary_lines.append("")
    summary_lines.append("--- Quality Metrics ---")
    for line in metric_sections["quality"]:
        summary_lines.append(line)
    summary_lines.append("")
    summary_lines.append("--- 16S Status ---")
    for line in metric_sections["rrna16s"]:
        summary_lines.append(line)
    summary_lines.append("")
    summary_lines.append("--- Redundancy ---")
    summary_lines.append(f"Exact duplicate CSV: {os.path.abspath(exact_duplicate_result['output_csv'])}")
    summary_lines.append(f"Exact duplicate groups: {exact_duplicate_result['duplicate_group_count']}")
    summary_lines.append(f"Genomes in exact duplicate groups: {exact_duplicate_result['duplicate_genome_count']}")
    if args.ani_cluster:
        summary_lines.append(f"ANI clustering status: {ani_cluster_result.get('status', 'unknown')}")
        if ani_cluster_result.get("message"):
            summary_lines.append(f"ANI clustering message: {ani_cluster_result['message']}")
        if ani_cluster_result.get("edge_csv"):
            summary_lines.append(f"ANI edge CSV: {os.path.abspath(str(ani_cluster_result['edge_csv']))}")
        if ani_cluster_result.get("cluster_csv"):
            summary_lines.append(f"ANI cluster CSV: {os.path.abspath(str(ani_cluster_result['cluster_csv']))}")
        if ani_cluster_result.get("log_path"):
            summary_lines.append(f"ANI log: {os.path.abspath(str(ani_cluster_result['log_path']))}")
        summary_lines.append(f"ANI cluster count: {ani_cluster_result.get('cluster_count', 0)}")
        summary_lines.append(f"Multi-member ANI clusters: {ani_cluster_result.get('multi_member_clusters', 0)}")
        summary_lines.append(f"Genomes in ANI multi-member clusters: {ani_cluster_result.get('clustered_genomes', 0)}")
        summary_lines.append(f"Qualifying ANI edges: {ani_cluster_result.get('qualifying_edges', 0)}")
    summary_lines.append("")
    summary_lines.append("--- Output ---")
    summary_lines.append(f"Standardized genomes: {len(records)}")
    summary_lines.append(f"Output dir: {os.path.abspath(base_out_dir)}")
    summary_lines.append(f"Genomes dir: {os.path.abspath(genomes_out_dir)}")
    summary_lines.append(f"Metadata CSV: {os.path.abspath(metadata_path)}")
    summary_lines.append(f"ID mapping CSV: {os.path.abspath(mapping_path)}")
    summary_lines.append(f"Scan inventory CSV: {os.path.abspath(scan_path)}")
    summary_lines.append(f"Source counts CSV: {os.path.abspath(source_inventory_path)}")
    summary_lines.append(f"Version comparison CSV: {os.path.abspath(version_compare_path)}")
    summary_lines.append("=" * 60)
    log_path = os.path.join(base_out_dir, "build_summary.log")
    write_summary_log(log_path, summary_lines)

    print("== Build Metrics ==")
    for line in metric_sections["size"]:
        print(f"  {line}")
    for line in metric_sections["quality"]:
        print(f"  {line}")
    if args.run_barrnap or args.collect_16s_to:
        for line in metric_sections["rrna16s"]:
            print(f"  {line}")
    print("== Redundancy ==")
    print(
        "  Exact duplicates: "
        f"{exact_duplicate_result['duplicate_group_count']} groups / "
        f"{exact_duplicate_result['duplicate_genome_count']} genomes"
    )
    if args.ani_cluster:
        print(
            "  ANI95 AF60 clusters: "
            f"status={ani_cluster_result.get('status')}, "
            f"clusters={ani_cluster_result.get('cluster_count', 0)}, "
            f"multi_member={ani_cluster_result.get('multi_member_clusters', 0)}"
        )

    if os.path.isdir(new_genomes_out_dir):
        shutil.rmtree(new_genomes_out_dir, ignore_errors=True)
    print(f"Done. Output: {os.path.abspath(base_out_dir)}")


def run(argv: Optional[List[str]] = None) -> int:
    """Entry point that allows both CLI invocation and import-time use.

    ``argv`` is reserved for future direct argv passthrough; the current
    implementation parses from ``sys.argv``.
    """
    del argv  # currently unused
    try:
        main()  # type: ignore[arg-type]
    except SystemExit as exc:
        return int(exc.code) if isinstance(exc.code, int) else 1
    return 0


if __name__ == "__main__":
    raise SystemExit(run())
