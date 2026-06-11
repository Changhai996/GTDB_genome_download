#!/usr/bin/env bash
# GTDB Renew — one-shot launcher.
#
# Usage:
#   ./run.sh build             --db-name Bathy --db-version v1 --database-root /path
#   ./run.sh download-gtdb     --release 220.0 --taxon "p__Bathyarchaeota"
#   ./run.sh run-all           --db-name Bathy --db-version v1 --release 220.0 \
#                              --taxon "p__Bathyarchaeota" --run-checkm2 --run-barrnap
#
# This script:
#   1. Ensures `pixi` is installed (downloads to ./.pixi/bin if missing).
#   2. Runs `pixi install` on first use.
#   3. Forwards every argument to `pixi run gttdb <args>`.
#
# Linux & macOS are supported. Tested on macOS 14 (arm64) and Linux x86_64.
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$PROJECT_ROOT"

# -- 1. Install pixi if absent ----------------------------------------------------
PIXI="$PROJECT_ROOT/.pixi/bin/pixi"
if [ ! -x "$PIXI" ]; then
  echo "[run.sh] pixi not found locally; downloading..."
  case "$(uname -s)" in
    Linux*)  TARGET_OS=linux ;;
    Darwin*) TARGET_OS=osx ;;
    *) echo "[run.sh] Unsupported OS: $(uname -s)" >&2; exit 1 ;;
  esac
  case "$(uname -m)" in
    x86_64|amd64)  TARGET_ARCH=x86_64 ;;
    aarch64|arm64) TARGET_ARCH=arm64 ;;
    *) echo "[run.sh] Unsupported arch: $(uname -m)" >&2; exit 1 ;;
  esac
  mkdir -p "$PROJECT_ROOT/.pixi/bin"
  curl -fsSL "https://pixi.sh/install.sh" | PIXI_HOME="$PROJECT_ROOT/.pixi" \
    PIXI_VERSION=latest bash -s -- -y
fi

# -- 2. Install dependencies on first run -----------------------------------------
if [ ! -d "$PROJECT_ROOT/.pixi/envs/default" ]; then
  echo "[run.sh] First run: resolving dependencies via pixi install ..."
  "$PIXI" install
fi

# -- 3. Forward to gttdb ----------------------------------------------------------
SUBCMD="${1:-}"
if [ -z "$SUBCMD" ]; then
  cat <<'USAGE'
GTDB Renew — launcher

Available commands:
  build            Build a versioned genome database from local source folders.
  download-gtdb    Download a GTDB release and filter genomes by taxon.
  run-all          Download a GTDB taxon and build a versioned database in one step.
  help             Show gttdb help.

Examples:
  ./run.sh build --db-name Bathy --db-version v1 --database-root /path/to/Database
  ./run.sh download-gtdb --release 220.0 --taxon "p__Bathyarchaeota"
  ./run.sh run-all --db-name Bathy --db-version v1 --release 220.0 \
       --taxon "p__Bathyarchaeota" --run-checkm2 --run-barrnap
USAGE
  exit 0
fi
shift || true

exec "$PIXI" run gttdb "$SUBCMD" "$@"
