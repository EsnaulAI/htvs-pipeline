#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SNAKEFILE="${ROOT_DIR}/workflow/Snakefile"

exec snakemake --unlock -s "${SNAKEFILE}" "$@"
