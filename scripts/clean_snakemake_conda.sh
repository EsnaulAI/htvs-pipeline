#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONDA_CACHE="${ROOT_DIR}/.snakemake/conda"

if [ -d "${CONDA_CACHE}" ]; then
  echo "Removing Snakemake conda cache: ${CONDA_CACHE}"
  rm -rf "${CONDA_CACHE}"
else
  echo "Snakemake conda cache not found at ${CONDA_CACHE}"
fi

echo "Done. Re-run Snakemake with --use-conda."
