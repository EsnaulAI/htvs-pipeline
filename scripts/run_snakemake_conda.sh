#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PROFILE_DIR="${ROOT_DIR}/workflow/profiles/conda"

if [ ! -d "${PROFILE_DIR}" ]; then
  echo "Profile not found at ${PROFILE_DIR}" >&2
  exit 1
fi

exec snakemake --profile "${PROFILE_DIR}" "$@"
