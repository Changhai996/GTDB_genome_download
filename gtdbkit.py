"""GTDB Renew — unified command-line interface.

Recommended subcommands
-----------------------

* ``web``         — open the Streamlit web interface.
* ``fetch``       — download GTDB genomes by taxon or accession.
* ``build-db``    — build/manage a local versioned genome database.
* ``prepare-db``  — convenience workflow: fetch GTDB genomes, then build DB.

Legacy aliases are kept for compatibility:
``build`` -> ``build-db``
``download-gtdb`` -> ``fetch``
``run-all`` -> ``prepare-db``
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys


def _add_build_parser(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser(
        "build-db",
        aliases=["build"],
        help="第三部分：构建和维护本地版本化数据库。",
        description=(
            "第三部分：扫描一个或多个数据源目录，统一重命名基因组与序列 header，"
            "提取 rRNA，按需运行 CheckM2，并输出版本化数据库与 metadata。"
            "原始文件不会被修改。"
        ),
    )
    p.add_argument("-n", "--name", "--db-name", dest="db_name", required=True, help="数据库名称。")
    p.add_argument("-v", "--version", "--db-version", dest="db_version", required=True, help="数据库版本号。")
    p.add_argument(
        "-D",
        "--source-root",
        "--database-root",
        dest="database_root",
        help="总输入目录。其下每个一级子目录视为一个来源。",
    )
    p.add_argument(
        "-s",
        "--source-dir",
        "--sources",
        dest="sources",
        nargs="+",
        help="一个或多个显式来源目录，可替代 --source-root。",
    )
    p.add_argument("-o", "--out-dir", "--output-root", dest="output_root", default="local_databases", help="输出目录。")
    p.add_argument("-j", "--threads", type=int, default=8, help="并行线程数。")
    p.add_argument("-Q", "--checkm2", "--run-checkm2", dest="run_checkm2", action="store_true", default=False, help="启用 CheckM2 质控。")
    p.add_argument("-B", "--barrnap", "--run-barrnap", dest="run_barrnap", action="store_true", default=False, help="启用 barrnap rRNA 提取。")
    p.add_argument(
        "-k",
        "--rrna-kingdom",
        "--barrnap-kingdom",
        dest="barrnap_kingdom",
        default="bac",
        choices=["bac", "arc", "euk", "mito"],
        help="barrnap kingdom 参数。",
    )
    p.add_argument(
        "-c",
        "--checkm2-db",
        "--checkm2-db-path",
        dest="checkm2_db_path",
        default=None,
        help="CheckM2 DIAMOND 数据库路径，可为 .dmnd 文件或其所在目录。",
    )
    p.add_argument(
        "--exclude-hidden",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="是否跳过隐藏文件，例如 .DS_Store 和 ._* AppleDouble 文件。默认开启。",
    )
    p.add_argument(
        "--strict-fasta-check",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="是否按文件内容校验 FASTA 格式，并自动识别无标准后缀的 FASTA。默认开启。",
    )
    p.add_argument(
        "-S",
        "--collect-16s-to",
        dest="collect_16s_to",
        default=None,
        help="把数据库里所有 16S 序列汇总到一个 FASTA 文件。",
    )
    p.add_argument("--ani-cluster", action=argparse.BooleanOptionalAction, default=False, help="启用 skani 的 ANI 重复度聚类分析。默认关闭。")
    p.add_argument("--ani-threshold", type=float, default=95.0, help="ANI 聚类阈值，默认 95.0。")
    p.add_argument("--af-threshold", type=float, default=60.0, help="AF 聚类阈值，默认 60.0。")


def _add_download_parser(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser(
        "fetch",
        aliases=["download-gtdb", "download"],
        help="第二部分：按 taxon 或 accession 下载 GTDB/NCBI 基因组。",
        description=(
            "第二部分：支持两种下载入口。1) 提供 GTDB taxon，从 GTDB metadata 解析 accession；"
            "2) 直接提供 accession。随后使用 NCBI datasets 分批并行下载。"
        ),
    )
    p.add_argument("-R", "--release", default="220.0", help="GTDB release，例如 220.0。")
    p.add_argument(
        "-t",
        "--taxon",
        action="append",
        help="GTDB taxon，可重复使用，也支持逗号分隔。",
    )
    p.add_argument("-a", "--accession", action="append", help="直接给 assembly accession，可重复或逗号分隔。")
    p.add_argument("-A", "--accession-file", default=None, help="包含 accession 列表的文本文件。")
    p.add_argument("-r", "--rank", default=None, choices=["d", "p", "c", "o", "f", "g", "s"], help="限定 GTDB rank。")
    p.add_argument(
        "-m",
        "--scope",
        "--mode",
        dest="mode",
        choices=["representative", "all", "accessions-only"],
        default="representative",
        help="下载范围：representative / all / accessions-only。",
    )
    p.add_argument("-o", "--out-dir", "--output-root", dest="output_root", default="gtdb_downloads", help="下载输出目录。")
    p.add_argument("-j", "--threads", type=int, default=4, help="并行下载任务数。")
    p.add_argument("-b", "--batch-size", type=int, default=50, help="每个 datasets 批次的 accession 数量。")
    p.add_argument(
        "-i",
        "--import-dir",
        "--import-to-dir",
        dest="import_to_dir",
        default=None,
        help="把下载得到的 fasta 直接导入指定目录，供第三部分使用。",
    )
    p.add_argument("--import-mode", choices=["symlink", "copy"], default="symlink", help="导入到 --import-dir 的方式。")
    p.add_argument("-L", "--list-releases", action="store_true", help="列出可用 GTDB release。")
    p.add_argument("--no-cache", action="store_true", help="不使用本地 GTDB metadata 缓存。")


def _add_prepare_parser(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser(
        "prepare-db",
        aliases=["run-all", "fetch-build"],
        help="第二+三部分：先下载，再导入并构建数据库。",
        description=(
            "先执行第二部分下载，再把结果与额外本地来源一起送入第三部分，"
            "完成数据库构建。这里需要 taxon/accession 只是因为第一步要先下载数据，"
            "并不是第三部分本身需要 taxon。"
        ),
    )
    p.add_argument("-n", "--name", "--db-name", dest="db_name", required=True, help="数据库名称。")
    p.add_argument("-v", "--version", "--db-version", dest="db_version", required=True, help="数据库版本号。")
    p.add_argument("-R", "--release", default="220.0", help="GTDB release。")
    p.add_argument("-t", "--taxon", action="append", help="要下载的 GTDB taxon。")
    p.add_argument("-a", "--accession", action="append", help="要下载的 accession。")
    p.add_argument("-A", "--accession-file", default=None, help="包含 accession 列表的文本文件。")
    p.add_argument("-r", "--rank", default=None, choices=["d", "p", "c", "o", "f", "g", "s"], help="限定 GTDB rank。")
    p.add_argument(
        "-m",
        "--scope",
        "--mode",
        dest="mode",
        choices=["representative", "all", "accessions-only"],
        default="representative",
        help="下载范围：representative / all / accessions-only。",
    )
    p.add_argument("-j", "--threads", type=int, default=8, help="下载和后续处理的并行线程数。")
    p.add_argument("-b", "--batch-size", type=int, default=50, help="每个下载批次的 accession 数量。")
    p.add_argument("-Q", "--checkm2", "--run-checkm2", dest="run_checkm2", action="store_true", default=False, help="启用 CheckM2。")
    p.add_argument("-B", "--barrnap", "--run-barrnap", dest="run_barrnap", action="store_true", default=False, help="启用 barrnap。")
    p.add_argument("-k", "--rrna-kingdom", "--barrnap-kingdom", dest="barrnap_kingdom", default="bac", choices=["bac", "arc", "euk", "mito"])
    p.add_argument("-c", "--checkm2-db", "--checkm2-db-path", dest="checkm2_db_path", default=None)
    p.add_argument(
        "--exclude-hidden",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="是否跳过隐藏文件，例如 .DS_Store 和 ._* AppleDouble 文件。默认开启。",
    )
    p.add_argument(
        "--strict-fasta-check",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="是否按文件内容校验 FASTA 格式，并自动识别无标准后缀的 FASTA。默认开启。",
    )
    p.add_argument(
        "-S",
        "--collect-16s-to",
        dest="collect_16s_to",
        default=None,
        help="把数据库里所有 16S 序列汇总到一个 FASTA 文件。",
    )
    p.add_argument("--ani-cluster", action=argparse.BooleanOptionalAction, default=False, help="启用 skani 的 ANI 重复度聚类分析。默认关闭。")
    p.add_argument("--ani-threshold", type=float, default=95.0, help="ANI 聚类阈值，默认 95.0。")
    p.add_argument("--af-threshold", type=float, default=60.0, help="AF 聚类阈值，默认 60.0。")
    p.add_argument(
        "-s",
        "--source-dir",
        "--extra-sources",
        dest="extra_sources",
        nargs="+",
        default=[],
        help="额外合并进数据库的本地来源目录。",
    )
    p.add_argument(
        "-d",
        "--download-dir",
        "--download-root",
        dest="download_root",
        default="gtdb_downloads",
        help="下载阶段输出目录。",
    )
    p.add_argument(
        "-o",
        "--out-dir",
        "--output-root",
        dest="output_root",
        default="local_databases",
        help="最终数据库输出目录。",
    )
    p.add_argument(
        "-i",
        "--import-dir",
        "--import-to-dir",
        dest="import_to_dir",
        default=None,
        help="把下载到的 GTDB fasta 导入指定目录，再作为第三部分输入。",
    )
    p.add_argument("--import-mode", choices=["symlink", "copy"], default="symlink")


def _add_web_parser(sub: argparse._SubParsersAction) -> None:
    sub.add_parser(
        "web",
        aliases=["ui"],
        help="第一部分：打开网页版 GTDB 比较与下载界面。",
        description="启动 Streamlit 网页界面。",
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


def _route_prepare(argv: list) -> int:
    # Parse just enough to call the two phases.
    parser = argparse.ArgumentParser(
        prog="gtdbkit prepare-db",
        description="先下载 GTDB 数据，再导入并构建本地版本化数据库。",
    )
    parser.add_argument("-n", "--name", "--db-name", dest="db_name", required=True)
    parser.add_argument("-v", "--version", "--db-version", dest="db_version", required=True)
    parser.add_argument("-R", "--release", default="220.0")
    parser.add_argument("-t", "--taxon", action="append")
    parser.add_argument("-a", "--accession", action="append")
    parser.add_argument("-A", "--accession-file", default=None)
    parser.add_argument("-r", "--rank", default=None, choices=["d", "p", "c", "o", "f", "g", "s"])
    parser.add_argument("-m", "--scope", "--mode", dest="mode", choices=["representative", "all", "accessions-only"], default="representative")
    parser.add_argument("-j", "--threads", type=int, default=8)
    parser.add_argument("-b", "--batch-size", type=int, default=50)
    parser.add_argument("-Q", "--checkm2", "--run-checkm2", dest="run_checkm2", action="store_true", default=False)
    parser.add_argument("-B", "--barrnap", "--run-barrnap", dest="run_barrnap", action="store_true", default=False)
    parser.add_argument("-k", "--rrna-kingdom", "--barrnap-kingdom", dest="barrnap_kingdom", default="bac", choices=["bac", "arc", "euk", "mito"])
    parser.add_argument("-c", "--checkm2-db", "--checkm2-db-path", dest="checkm2_db_path", default=None)
    parser.add_argument("--exclude-hidden", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--strict-fasta-check", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("-S", "--collect-16s-to", dest="collect_16s_to", default=None)
    parser.add_argument("--ani-cluster", action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument("--ani-threshold", type=float, default=95.0)
    parser.add_argument("--af-threshold", type=float, default=60.0)
    parser.add_argument("-s", "--source-dir", "--extra-sources", dest="extra_sources", nargs="+", default=[])
    parser.add_argument("-d", "--download-dir", "--download-root", dest="download_root", default="gtdb_downloads")
    parser.add_argument("-o", "--out-dir", "--output-root", dest="output_root", default="local_databases")
    parser.add_argument("-i", "--import-dir", "--import-to-dir", dest="import_to_dir", default=None)
    parser.add_argument("--import-mode", choices=["symlink", "copy"], default="symlink")
    args = parser.parse_args(argv)

    if not args.taxon and not args.accession and not args.accession_file:
        raise SystemExit("prepare-db requires --taxon or --accession/--accession-file.")

    dl_argv = ["--release", args.release, "--mode", args.mode,
               "--output-root", args.download_root, "--threads", str(args.threads),
               "--batch-size", str(args.batch_size)]
    if args.taxon:
        dl_argv += ["--taxon", *args.taxon]
    if args.accession:
        dl_argv += ["--accession", *args.accession]
    if args.accession_file:
        dl_argv += ["--accession-file", args.accession_file]
    if args.rank:
        dl_argv += ["--rank", args.rank]
    if args.import_to_dir:
        dl_argv += ["--import-to-dir", args.import_to_dir, "--import-mode", args.import_mode]
    print("== Phase 1: download GTDB ==", flush=True)
    rc = _route_download(dl_argv)
    if rc != 0:
        return rc

    # Compose sources: every per-taxon fasta dir we just created.
    sources = list(args.extra_sources or [])
    if args.import_to_dir and os.path.isdir(args.import_to_dir):
        sources.append(args.import_to_dir)
    by_taxon = os.path.join(args.download_root, "by_taxon")
    if os.path.isdir(by_taxon):
        for taxon_dir in sorted(os.listdir(by_taxon)):
            fasta_dir = os.path.join(by_taxon, taxon_dir, "fasta")
            if os.path.isdir(fasta_dir) and any(f.endswith((".fna", ".fa", ".fasta", ".fna.gz", ".fa.gz", ".fasta.gz")) for f in os.listdir(fasta_dir)):
                sources.append(fasta_dir)
    by_accession = os.path.join(args.download_root, "by_accession")
    if os.path.isdir(by_accession):
        for query_dir in sorted(os.listdir(by_accession)):
            fasta_dir = os.path.join(by_accession, query_dir, "fasta")
            if os.path.isdir(fasta_dir) and any(f.endswith((".fna", ".fa", ".fasta", ".fna.gz", ".fa.gz", ".fasta.gz")) for f in os.listdir(fasta_dir)):
                sources.append(fasta_dir)
    if not sources:
        print("prepare-db: no source folders available after download; aborting build phase.", flush=True)
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
    if not args.exclude_hidden:
        build_argv += ["--no-exclude-hidden"]
    if not args.strict_fasta_check:
        build_argv += ["--no-strict-fasta-check"]
    if args.collect_16s_to:
        build_argv += ["--collect-16s-to", args.collect_16s_to]
    if args.ani_cluster:
        build_argv += ["--ani-cluster"]
    if args.ani_threshold != 95.0:
        build_argv += ["--ani-threshold", str(args.ani_threshold)]
    if args.af_threshold != 60.0:
        build_argv += ["--af-threshold", str(args.af_threshold)]

    print("\n== Phase 2: build versioned database ==", flush=True)
    rc = _route_build(build_argv)
    return rc


def _route_web(argv: list) -> int:
    app_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "app.py")
    cmd = [sys.executable, "-m", "streamlit", "run", app_path] + list(argv)
    completed = subprocess.run(cmd)
    return completed.returncode


def main() -> int:
    p = argparse.ArgumentParser(
        prog="gtdbkit",
        description=(
            "GTDB Renew 统一入口：web = 网页版；fetch = GTDB 下载；"
            "build-db = 本地数据库构建；prepare-db = 下载后直接入库。"
        ),
    )
    sub = p.add_subparsers(dest="cmd", required=True)
    _add_build_parser(sub)
    _add_download_parser(sub)
    _add_prepare_parser(sub)
    _add_web_parser(sub)
    if len(sys.argv) < 2 or sys.argv[1] in {"-h", "--help"}:
        p.print_help()
        return 0
    cmd = sys.argv[1]
    rest = sys.argv[2:]
    if cmd in {"build-db", "build"}:
        return _route_build(rest)
    if cmd in {"fetch", "download-gtdb", "download"}:
        return _route_download(rest)
    if cmd in {"prepare-db", "run-all", "fetch-build"}:
        return _route_prepare(rest)
    if cmd in {"web", "ui"}:
        return _route_web(rest)
    p.error(f"Unknown subcommand: {cmd}")


if __name__ == "__main__":
    raise SystemExit(main())
