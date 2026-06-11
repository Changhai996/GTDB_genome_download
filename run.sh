#!/usr/bin/env bash
# GTDB Renew — one-shot launcher.
#
# Usage:
#   ./run.sh web
#   ./run.sh fetch            -R 220.0 -t "p__Bathyarchaeota"
#   ./run.sh build-db         -n Bathy -v v1 -D /path
#   ./run.sh prepare-db       -n Bathy -v v1 -t "p__Bathyarchaeota" -Q -B
#
# This script:
#   1. Ensures `pixi` is installed (downloads to ./.pixi/bin if missing).
#   2. Runs `pixi install` on first use.
#   3. Forwards every argument to `pixi run gtdbkit <args>`.
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

# -- 3. Forward to gtdbkit --------------------------------------------------------
SUBCMD="${1:-}"
if [ -z "$SUBCMD" ]; then
  cat <<'USAGE'
GTDB Renew — launcher

Available commands:
  web              Open the Streamlit web UI.
  fetch            Download GTDB genomes by taxon or accession.
  build-db         Build a versioned genome database from local source folders.
  prepare-db       Download GTDB genomes and then build a database.
  help             Show gtdbkit help.

Examples:
  ./run.sh web
  ./run.sh fetch -R 220.0 -t "p__Bathyarchaeota"
  ./run.sh build-db -n Bathy --version v1 -D /path/to/Database
  ./run.sh prepare-db -n Bathy -v v1 -t "p__Bathyarchaeota" -Q -B
USAGE
  exit 0
fi
shift || true

if [ "$SUBCMD" = "help" ]; then
  exec "$PIXI" run gtdbkit --help
fi

exec "$PIXI" run gtdbkit "$SUBCMD" "$@"
