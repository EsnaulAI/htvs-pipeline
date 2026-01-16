#!/usr/bin/env bash
set -euo pipefail

IMAGE="${GNINA_CONTAINER_IMAGE:-gnina/gnina}"
RUNTIME="${GNINA_CONTAINER_RUNTIME:-}"
HOST_DIR="${GNINA_CONTAINER_HOST_DIR:-$PWD}"
WORKDIR="${GNINA_CONTAINER_WORKDIR:-/work}"

if [[ -z "$RUNTIME" ]]; then
  if command -v docker >/dev/null 2>&1; then
    RUNTIME="docker"
  elif command -v apptainer >/dev/null 2>&1; then
    RUNTIME="apptainer"
  elif command -v singularity >/dev/null 2>&1; then
    RUNTIME="singularity"
  else
    echo "No se encontró Docker, Apptainer o Singularity para ejecutar Gnina." >&2
    exit 127
  fi
fi

case "$RUNTIME" in
  docker)
    exec docker run --rm \
      -v "${HOST_DIR}:${WORKDIR}" \
      -w "${WORKDIR}" \
      "$IMAGE" \
      gnina "$@"
    ;;
  apptainer|singularity)
    exec "$RUNTIME" exec \
      --bind "${HOST_DIR}:${WORKDIR}" \
      "$IMAGE" \
      gnina "$@"
    ;;
  *)
    echo "Runtime no soportado: ${RUNTIME}" >&2
    exit 2
    ;;
esac
