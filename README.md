# GTDB Renew

End-to-end pipeline for:

1. **Downloading GTDB genomes** for a specified release, filtered to any
   GTDB taxon (phylum, class, order, family, genus, species). The downloader
   reads the GTDB taxonomy + metadata tables, extracts assembly accession
   numbers, and then launches multiple parallel `datasets download genome accession`
   tasks. No full GTDB FASTA tar package is downloaded.
2. **Building a versioned local genome database** from those genomes
   (and any additional local sources). Each version is reproducible:
   standardized filenames, contig-header renaming to match the genome
   name, basic statistics, optional barrnap rRNA extraction, optional
   CheckM2 quality assessment.

The project is fully managed by [Pixi](https://pixi.sh) — one command
resolves Python 3.12, `checkm2`, `barrnap`, `pandas`, `requests`, and
the DIAMOND runtime.

## Supported platforms

| Platform | Status |
|---|---|
| macOS 14 (arm64) | ✅ tested |
| Linux x86_64     | ✅ tested |
| Linux aarch64    | ⏳ should work (bioconda checkm2 available) |
| Windows          | ❌ not supported (checkm2 + multiprocessing issues) |

## One-shot install (Linux)

```bash
git clone <your-fork-url> GTDB_renew
cd GTDB_renew
chmod +x run.sh
./run.sh run-all \
    --db-name Bathyarchaeia \
    --db-version v1.0 \
    --release 220.0 \
    --taxon "p__Bathyarchaeota" \
    --run-checkm2 \
    --run-barrnap
```

The first run downloads Pixi (~10 MB) and resolves ~3 GB of
dependencies (CheckM2 + DIAMOND + barrnap + Python 3.12 + friends).
Subsequent runs are fast — Pixi caches everything.

## Subcommands

### `download-gtdb`

```bash
./run.sh download-gtdb --release 220.0 --taxon "p__Bathyarchaeota" --mode representative
./run.sh download-gtdb --release 220.0 --taxon "c__Bathyarchaeia" --mode all
./run.sh download-gtdb --list-releases
```

- `--mode representative` keeps only rows where GTDB metadata marks the
  genome as the GTDB representative, extracts the assembly accession from
  the metadata table, then downloads those genomes with batched
  `datasets download genome accession` jobs.
- `--mode all` keeps every matched GTDB genome, extracts their assembly
  accession numbers from metadata, and downloads them the same way.
- `--mode accessions-only` stops after writing `accessions.csv`, which is
  useful for auditing the selected accessions before any download starts.
- `--taxon` is repeatable and accepts any GTDB rank
  (`p__`, `c__`, `o__`, `f__`, `g__`, `s__`).
- `--rank {d|p|c|o|f|g|s}` is optional and forces the rank of
  matching taxa.
- `--batch-size` controls how many accessions are bundled into each
  `datasets` task; `--threads` controls how many batches run concurrently.

### `build`

```bash
./run.sh build \
    --db-name Bathyarchaeia \
    --db-version v1.0 \
    --database-root /path/to/Database \
    --output-root ./local_databases \
    --threads 8 \
    --run-checkm2 \
    --run-barrnap \
    --checkm2-db-path /path/to/CheckM2_database/uniref100.KO.1.dmnd
```

Or specify sources directly:

```bash
./run.sh build \
    --db-name Bathyarchaeia --db-version v1.0 \
    --sources /path/to/Source_A /path/to/Source_B \
    --threads 8 --run-checkm2 --run-barrnap
```

Each invocation:

- Walks `--database-root/<source>/...` (or every `--sources`) and
  resolves FASTA files at any depth.
- Copies each FASTA into
  `<output-root>/<db-name>_<db-version>/genomes/<standardized>.fna`
  with **contig headers renamed to match the standardized genome name**.
  Original files are never modified.
- Records `Original_Name / Original_Path / Source_Folder /
  Standardized_Name / File_MD5` in a mapping CSV.
- Computes `Genome_Size_bp / Contig_Count / GC_Content_%`.
- Runs **barrnap** to extract rRNA and writes
  `<standardized>.gff`, `<standardized>_rrna.fasta` (the union of all
  rRNA types), and `<standardized>_16S.fasta` (16S slices taken
  directly from the genome using GFF coordinates — independent of
  barrnap's own `--outseq`).
- Runs **CheckM2** on the new genomes, reports
  `Completeness / Contamination / GTDB_Score`, and folds
  `Quality_Status` (`ok` / `no_database` / `failed` / `skipped` /
  `carried_forward`) into the metadata.
- Compares to the previous version snapshot and marks new vs.
  carried-forward genomes.
- Writes a comprehensive `build_summary.log` with parameter
  snapshot, source counts, CheckM2 status, and file inventory.

### `run-all`

Convenience wrapper: `download-gtdb` followed by `build`.

```bash
./run.sh run-all \
    --db-name Bathyarchaeia \
    --db-version v1.0 \
    --release 220.0 \
    --taxon "p__Bathyarchaeota" \
    --run-checkm2 \
    --run-barrnap \
    --checkm2-db-path /path/to/CheckM2_database/uniref100.KO.1.dmnd
```

Use `--extra-sources /path/A /path/B` to merge local genomes with the
GTDB pull.

## Incremental versioning

The build step looks for a previous snapshot of the same `db-name` and
copies the existing `genomes/` + `barrnap_results/` into the new
output directory, only processing **new** files (MD5-based). To force
a re-build, change `--db-version` (e.g. `v1.0` → `v1.1`).

## CheckM2 caveat (Python 3.12)

CheckM2 1.1.0 has a known multiprocessing pickling bug under
Python 3.12 (`AttributeError: 'Predictor' object has no attribute
'__set_up_prodigal_thread'`). This project ships two complementary
fixes:

1. `checkm2_compat/sitecustomize.py` — auto-injected via
   `PYTHONPATH` so that `multiprocessing.Pool` uses the `fork`
   start method.
2. An in-place patch at the bottom of
   `<pixi-env>/lib/python3.12/site-packages/checkm2/predictQuality.py`
   that adds unmangled aliases for the four private `Predictor`
   methods. This patch is applied automatically by `db_builder_cli.py`
   the first time `--run-checkm2` succeeds — look for the comment
   "GTDB Renew compatibility patch" at the bottom of the file. If
   you `pixi install` again and the patch is lost, simply re-run
   `./run.sh build --run-checkm2` once; on failure the patch is
   reapplied.

If you would rather avoid touching `site-packages`, set
`--threads 1` and the bug disappears — slower, but safe.

## CheckM2 database path

The build step resolves the DIAMOND database in this order:

1. `--checkm2-db-path` (a `.dmnd` file or a directory containing one)
2. `$CHECKM2DB` / `$DIAMOND_DB` environment variable
3. `<project>/checkm2_database/CheckM2_database/*.dmnd`
4. `<project>/checkm2_database/*.dmnd`

The first one that exists wins. The resolved path is written into the
CheckM2 child-process environment automatically.

Download a CheckM2 database with:

```bash
pixi run checkm2 database --download
```

The default download location is `~/databases/CheckM2_database/`;
move it into `<project>/checkm2_database/CheckM2_database/` to use
the project-local default lookup, or pass it via
`--checkm2-db-path`.

## Outputs of a single build

```
local_databases/Bathyarchaeia_v1.0/
├── genomes/                          standardized per-genome fasta
├── barrnap_results/                  *.gff, *_rrna.fasta, *_16S.fasta
├── checkm2_results/                  quality_report.tsv, checkm2.log
├── Bathyarchaeia_v1.0_metadata.csv   one row per genome
├── Bathyarchaeia_v1.0_genome_id_mapping.csv
├── Bathyarchaeia_v1.0_scan_inventory.csv
├── Bathyarchaeia_v1.0_source_counts.csv
├── Bathyarchaeia_v1.0_version_comparison.csv
└── build_summary.log
```

## Layout

```
GTDB_renew/
├── pixi.toml                        pixi workspace + tasks
├── gttdb.py                         unified CLI (subcommands)
├── db_builder_cli.py                versioned database builder
├── gtdb_downloader.py               GTDB downloader + taxon filter
├── checkm2_compat/
│   └── sitecustomize.py             Python 3.12 compat shim
├── run.sh                           one-shot launcher
└── README.md                        this file
```

## License

Internal tool. CheckM2, barrnap, and DIAMOND are governed by their
own licenses.
