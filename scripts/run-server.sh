#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
plugin_dir="$(cd -- "$script_dir/.." && pwd)"

if ! command -v Rscript >/dev/null 2>&1; then
  echo "r2-depmap requires Rscript on PATH. Run scripts/install-dependencies.R after installing R." >&2
  exit 127
fi

export R2_DEPMAP_PLUGIN_DIR="$plugin_dir"
export R2_DEPMAP_CONFIG="${R2_DEPMAP_CONFIG:-$plugin_dir/config/defaults.json}"
exec Rscript --vanilla "$plugin_dir/R/server.R"
