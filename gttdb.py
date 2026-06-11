"""GTDB Renew — unified command-line interface.

Subcommands
-----------

* ``build``            — build a versioned genome database from one or more
                         source folders (delegates to ``db_builder_cli``).
* ``download-gtdb``    — download a GTDB release and filter genomes by
                         taxon (delegates to ``gtdb_downloader``).
* ``run-all``          — convenience: download → build in a single step.

Use ``python -m gttdb <subcommand> --help`` (or ``pixi run gt <subcommand>
--help``) to see options.
"""
from __future__ import annotations

import argparse
import os
import sys


def _add_build_parser(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser(
        "build",
        help="Build a versioned genome database from local source folders.",
        description=(
            "Walks the source folders, renames genomes, extracts rRNA via "
            "barrnap, optionally runs CheckM2, and writes a versioned "
            "metadata bundle. Original files are never modified."
        ),
    )
    p.add_argument("--db-name", required=True)
    p.add_argument("--db-version", required=True)
    p.add_argument(
        "--database-root",
        help="Root directory containing one subfolder per source. Each top-level subfolder is treated as a separate source.",
    )
    p.add_argument(
        "--sources",
        nargs="+",
        help="One or more explicit source directories (alternative to --database-root).",
    )
    p.add_argument("--output-root", default="local_databases")
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--run-checkm2", action="store_true", default=False)
    p.add_argument("--run-barrnap", action="store_true", default=False)
    p.add_argument(
        "--barrnap-kingdom",
        default="bac",
        choices=["bac", "arc", "euk", "mito"],
    )
    p.add_argument(
        "--checkm2-db-path",
        default=None,
        help="Explicit path to the CheckM2 DIAMOND database file (*.dmnd) or to a directory containing one.",
    )


def _add_download_parser(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser(
        "download-gtdb",
        help="Download a GTDB release and filter genomes by taxon.",
        description=(
            "Downloads the GTDB taxonomy TSVs for the chosen release, "
            "filters by any GTDB taxon (phylum/class/order/family/genus/species), "
            "extracts assembly accession numbers from GTDB metadata, and "
            "launches multiple batched NCBI `datasets` download tasks in parallel."
        ),
    )
    p.add_argument("--release", default="220.0")
    p.add_argument(
        "--taxon",
        action="append",
        required=True,
        help="GTDB taxon to keep. Repeatable. Comma-separated values inside one --taxon are also accepted.",
    )
    p.add_argument("--rank", default=None, choices=["d", "p", "c", "o", "f", "g", "s"])
    p.add_argument(
        "--mode",
        choices=["representative", "all", "accessions-only"],
        default="representative",
    )
    p.add_argument("--output-root", default="gtdb_downloads")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--batch-size", type=int, default=50)
    p.add_argument("--list-releases", action="store_true")
    p.add_argument("--no-cache", action="store_true")


def _add_runall_parser(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser(
        "run-all",
        help="Download a GTDB taxon and build a versioned database in one step.",
        description=(
            "Convenience wrapper: 1) download GTDB genomes for the requested "
            "taxa, 2) build a versioned database from those genomes (and any "
            "additional --sources folders you pass)."
        ),
    )
    p.add_argument("--db-name", required=True)
    p.add_argument("--db-version", required=True)
    p.add_argument("--release", default="220.0")
    p.add_argument("--taxon", action="append", required=True)
    p.add_argument("--rank", default=None, choices=["d", "p", "c", "o", "f", "g", "s"])
    p.add_argument(
        "--mode",
        choices=["representative", "all", "accessions-only"],
        default="representative",
    )
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--batch-size", type=int, default=50)
    p.add_argument("--run-checkm2", action="store_true", default=False)
    p.add_argument("--run-barrnap", action="store_true", default=False)
    p.add_argument("--barrnap-kingdom", default="bac", choices=["bac", "arc", "euk", "mito"])
    p.add_argument("--checkm2-db-path", default=None)
    p.add_argument(
        "--extra-sources",
        nargs="+",
        default=[],
        help="Additional local source folders to merge with the GTDB download.",
    )
    p.add_argument(
        "--download-root",
        default="gtdb_downloads",
        help="Where to put intermediate GTDB downloads (per-taxon subfolders).",
    )
    p.add_argument(
        "--output-root",
        default="local_databases",
        help="Where to put the final versioned database snapshot.",
    )


def _route_build(argv: list) -> int:
    from db_builder_cli import main as build_main
    sys.argv = ["db_builder_cli"] + argv
    try:
        build_main()
    except SystemExit as exc:
        return int(exc.code) if isinstance(exc.code, int) else 1
    return 0


def _route_download(argv: list) -> int:
    from gtdb_downloader import main as dl_main
    sys.argv = ["gtdb_downloader"] + argv
    try:
        dl_main()
    except SystemExit as exc:
        return int(exc.code) if isinstance(exc.code, int) else 1
    return 0


def _route_runall(argv: list) -> int:
    # Parse just enough to call the two phases.
    parser = argparse.ArgumentParser(prog="gttdb run-all", add_help=False)
    parser.add_argument("--db-name", required=True)
    parser.add_argument("--db-version", required=True)
    parser.add_argument("--release", default="220.0")
    parser.add_argument("--taxon", action="append", required=True)
    parser.add_argument("--rank", default=None, choices=["d", "p", "c", "o", "f", "g", "s"])
    parser.add_argument("--mode", choices=["representative", "all", "accessions-only"], default="representative")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--batch-size", type=int, default=50)
    parser.add_argument("--run-checkm2", action="store_true", default=False)
    parser.add_argument("--run-barrnap", action="store_true", default=False)
    parser.add_argument("--barrnap-kingdom", default="bac", choices=["bac", "arc", "euk", "mito"])
    parser.add_argument("--checkm2-db-path", default=None)
    parser.add_argument("--extra-sources", nargs="+", default=[])
    parser.add_argument("--download-root", default="gtdb_downloads")
    parser.add_argument("--output-root", default="local_databases")
    args = parser.parse_args(argv)

    dl_argv = ["--release", args.release, "--taxon", *args.taxon, "--mode", args.mode,
               "--output-root", args.download_root, "--threads", str(args.threads),
               "--batch-size", str(args.batch_size)]
    if args.rank:
        dl_argv += ["--rank", args.rank]
    print("== Phase 1: download GTDB ==", flush=True)
    rc = _route_download(dl_argv)
    if rc != 0:
        return rc

    # Compose sources: every per-taxon fasta dir we just created.
    sources = list(args.extra_sources or [])
    by_taxon = os.path.join(args.download_root, "by_taxon")
    if os.path.isdir(by_taxon):
        for taxon_dir in sorted(os.listdir(by_taxon)):
            fasta_dir = os.path.join(by_taxon, taxon_dir, "fasta")
            if os.path.isdir(fasta_dir) and any(f.endswith((".fna", ".fa", ".fasta", ".fna.gz", ".fa.gz", ".fasta.gz")) for f in os.listdir(fasta_dir)):
                sources.append(fasta_dir)
    if not sources:
        print("run-all: no source folders available after download; aborting build phase.", flush=True)
        return 2

    # Database root = a temp dir that contains one subfolder per source.
    # Easier: pass --sources directly to the build phase.
    build_argv = [
        "--db-name", args.db_name,
        "--db-version", args.db_version,
        "--sources", *sources,
        "--output-root", args.output_root,
        "--threads", str(args.threads),
    ]
    if args.run_checkm2:
        build_argv += ["--run-checkm2"]
    if args.run_barrnap:
        build_argv += ["--run-barrnap", "--barrnap-kingdom", args.barrnap_kingdom]
    if args.checkm2_db_path:
        build_argv += ["--checkm2-db-path", args.checkm2_db_path]

    print("\n== Phase 2: build versioned database ==", flush=True)
    rc = _route_build(build_argv)
    return rc


def main() -> int:
    p = argparse.ArgumentParser(
        prog="gttdb",
        description="GTDB Renew: download GTDB genomes and build a versioned local database.",
    )
    sub = p.add_subparsers(dest="cmd", required=True)
    _add_build_parser(sub)
    _add_download_parser(sub)
    _add_runall_parser(sub)
    # Important: do NOT use parse_known_args() here. With subparsers, the
    # subparser's own parser consumes --db-name etc. and parse_known_args
    # on the root parser then returns an empty rest. We only need the
    # subcommand identifier from root and then take the raw remaining args
    # straight from sys.argv.
    if len(sys.argv) < 2:
        p.error("missing subcommand")
    cmd = sys.argv[1]
    rest = sys.argv[2:]
    if cmd == "build":
        return _route_build(rest)
    if cmd == "download-gtdb":
        return _route_download(rest)
    if cmd == "run-all":
        return _route_runall(rest)
    p.error(f"Unknown subcommand: {cmd}")


if __name__ == "__main__":
    raise SystemExit(main())
