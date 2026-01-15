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
        input=SimpleNamespace(
            receptor="results/module1/target_prepared.pdbqt",
            top_candidates="results/module3/top_candidates.csv",
        ),
        output=SimpleNamespace(scores="results/module3/gnina_scores.csv"),
        params=SimpleNamespace(pose_dir="results/module3/docking"),
        wildcards=SimpleNamespace(),
    )

import os
import re
import subprocess
from typing import Optional, Tuple

import pandas as pd

from logging_utils import (
    confirm_file,
    ensure_parent_dir,
    log_info,
    log_warn,
    require_file,
)


AFFINITY_RE = re.compile(r"Affinity:\s*(-?\d+(?:\.\d+)?)")
CNN_SCORE_RE = re.compile(r"CNNscore:\s*(-?\d+(?:\.\d+)?)", re.IGNORECASE)


def parse_gnina_output(text: str) -> Tuple[Optional[float], Optional[float]]:
    affinity = None
    cnn_score = None

    affinity_match = AFFINITY_RE.search(text)
    if affinity_match:
        affinity = float(affinity_match.group(1))

    cnn_match = CNN_SCORE_RE.search(text)
    if cnn_match:
        cnn_score = float(cnn_match.group(1))

    return affinity, cnn_score


def run_gnina(receptor: str, ligand_path: str) -> Tuple[Optional[float], Optional[float]]:
    command = [
        "gnina",
        "--receptor",
        receptor,
        "--ligand",
        ligand_path,
        "--score_only",
    ]
    result = subprocess.run(
        command,
        check=False,
        capture_output=True,
        text=True,
    )
    output = (result.stdout or "") + "\n" + (result.stderr or "")

    if result.returncode != 0:
        log_warn(
            f"⚠️ Gnina falló para {ligand_path} (código {result.returncode})."
        )
        return None, None

    affinity, cnn_score = parse_gnina_output(output)
    if affinity is None or cnn_score is None:
        log_warn(f"⚠️ No se detectaron scores completos en {ligand_path}.")

    return affinity, cnn_score


def load_ligand_ids(top_candidates_path: str) -> list[str]:
    df = pd.read_csv(top_candidates_path)
    if df.empty:
        return []
    if "ligand_id" in df.columns:
        return df["ligand_id"].dropna().astype(str).tolist()
    return df.iloc[:, 0].dropna().astype(str).tolist()


def main():
    receptor = snakemake.input.receptor
    top_candidates_path = snakemake.input.top_candidates
    output_csv = snakemake.output.scores
    pose_dir = getattr(snakemake.params, "pose_dir", "results/module3/docking")

    require_file(receptor, "receptor preparado")
    require_file(top_candidates_path, "candidatos top")

    ligand_ids = load_ligand_ids(top_candidates_path)
    ensure_parent_dir(output_csv)

    if not ligand_ids:
        log_warn("⚠️ No hay ligandos para re-score con Gnina.")
        pd.DataFrame(columns=["ligand_id", "vina_score", "cnn_score"]).to_csv(
            output_csv, index=False
        )
        confirm_file(output_csv, "scores Gnina")
        return

    records = []
    for ligand_id in ligand_ids:
        pose_path = os.path.join(pose_dir, f"{ligand_id}_out.pdbqt")
        if not os.path.isfile(pose_path):
            log_warn(f"⚠️ Pose no encontrada para {ligand_id}: {pose_path}")
            continue
        vina_score, cnn_score = run_gnina(receptor, pose_path)
        records.append(
            {
                "ligand_id": ligand_id,
                "vina_score": vina_score,
                "cnn_score": cnn_score,
            }
        )

    df = pd.DataFrame(records)
    df.to_csv(output_csv, index=False)
    confirm_file(output_csv, "scores Gnina")
    log_info(f"✅ Tabla Gnina lista: {output_csv}")


if __name__ == "__main__":
    main()
