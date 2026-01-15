# --- MOCK PARA DESARROLLO (Pylance no se quejará) ---
if "snakemake" not in globals():
    from pathlib import Path
    from types import SimpleNamespace
    import yaml

    def load_config():
        config_path = Path(__file__).resolve().parents[2] / "config" / "config.yaml"
        with open(config_path, "r") as f:
            return yaml.safe_load(f)

    config = load_config()
    snakemake = SimpleNamespace(
        input=SimpleNamespace(logs_dir="results/module3/logs"),
        output=SimpleNamespace(scores="results/module3/docking_scores.csv"),
        params=SimpleNamespace(pose_dir="results/module3/docking"),
        wildcards=SimpleNamespace(),
    )

import os
import re
import glob
import pandas as pd

from logging_utils import (
    confirm_file,
    ensure_parent_dir,
    log_info,
    log_warn,
)


AFFINITY_TABLE_RE = re.compile(r"^\s*(\d+)\s+(-?\d+(?:\.\d+)?)\s+")
AFFINITY_LINE_RE = re.compile(r"Affinity:\s*(-?\d+(?:\.\d+)?)")


def collect_log_files(inputs):
    items = list(inputs)
    log_files = []

    for item in items:
        if os.path.isdir(item):
            patterns = ["*.log", "*.txt", "*.out"]
            matched = []
            for pattern in patterns:
                matched.extend(glob.glob(os.path.join(item, pattern)))
            if not matched:
                matched = [
                    os.path.join(item, name)
                    for name in os.listdir(item)
                    if os.path.isfile(os.path.join(item, name))
                ]
            log_files.extend(matched)
        else:
            log_files.append(item)

    return sorted({path for path in log_files if os.path.isfile(path)})


def parse_best_affinity(log_path):
    scores = []
    with open(log_path, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            match = AFFINITY_TABLE_RE.match(line)
            if match:
                scores.append(float(match.group(2)))
                continue
            match = AFFINITY_LINE_RE.search(line)
            if match:
                scores.append(float(match.group(1)))

    if not scores:
        return None
    return min(scores)


def main():
    output_csv = getattr(snakemake.output, "scores", None) or snakemake.output[0]
    logs_input = getattr(snakemake.input, "logs", None)
    logs_dir = getattr(snakemake.input, "logs_dir", None)
    if logs_input is not None:
        input_items = list(logs_input)
    elif logs_dir is not None:
        input_items = [logs_dir]
    else:
        input_items = list(snakemake.input)
    log_paths = collect_log_files(input_items)
    pose_dir = getattr(snakemake.params, "pose_dir", "results/module3/docking")

    if not log_paths:
        log_warn("⚠️ No se encontraron logs de Smina para analizar.")

    ensure_parent_dir(output_csv)
    log_info(f"📥 Procesando {len(log_paths)} logs de Smina...")

    records = []
    for log_path in log_paths:
        affinity = parse_best_affinity(log_path)
        ligand_id = os.path.splitext(os.path.basename(log_path))[0]
        pose_path = os.path.join(pose_dir, f"{ligand_id}_out.pdbqt")
        if affinity is None:
            log_warn(f"⚠️ Sin afinidad detectada en: {log_path}")
            continue
        if not os.path.isfile(pose_path):
            log_warn(f"⚠️ No se encontró pose para {ligand_id}: {pose_path}")
        records.append(
            {
                "ligand_id": ligand_id,
                "score": affinity,
                "log_path": log_path,
                "pose_path": pose_path,
            }
        )

    df = pd.DataFrame(records)
    if not df.empty:
        df = df.sort_values("score", ascending=True)
        df = df[["ligand_id", "score", "log_path", "pose_path"]]
    else:
        df = pd.DataFrame(columns=["ligand_id", "score", "log_path", "pose_path"])

    df.to_csv(output_csv, index=False)
    confirm_file(output_csv, "scores de docking")
    log_info(f"✅ Tabla de scores generada: {output_csv}")


if __name__ == "__main__":
    main()
