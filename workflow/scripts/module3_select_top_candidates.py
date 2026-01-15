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
            scores="results/module3/docking_scores.csv",
            library="results/module2/library_enumerated.smi",
        ),
        output=SimpleNamespace(
            csv="results/module3/top_candidates.csv",
            smi="results/module3/top_candidates.smi",
        ),
        params=SimpleNamespace(
            top_percent=config.get("docking", {}).get("top_candidates_percent", 5)
        ),
        wildcards=SimpleNamespace(),
    )

import math
from typing import Dict, Tuple

import pandas as pd

from logging_utils import (
    confirm_file,
    ensure_parent_dir,
    log_info,
    log_warn,
    require_file,
)


def safe_ligand_id(name: str) -> str:
    safe_name = "".join([c for c in name if c.isalnum() or c in ("_", "-")])
    return safe_name or name


def load_library_smi(smi_path: str) -> Dict[str, Tuple[str, str]]:
    mapping: Dict[str, Tuple[str, str]] = {}
    with open(smi_path, "r", encoding="utf-8") as handle:
        for idx, line in enumerate(handle):
            if not line.strip():
                continue
            parts = line.split()
            smiles = parts[0]
            ligand_id = parts[1] if len(parts) > 1 else f"Ligand_{idx}"
            safe_id = safe_ligand_id(ligand_id)
            if safe_id in mapping:
                log_warn(f"Duplicado en SMILES: {safe_id} (se conserva el primero)")
                continue
            mapping[safe_id] = (smiles, ligand_id)
    return mapping


def clamp_top_percent(value: float) -> float:
    if value is None:
        return 5.0
    try:
        value = float(value)
    except (TypeError, ValueError):
        return 5.0
    return min(5.0, max(1.0, value))


def main():
    scores_path = snakemake.input.scores
    library_path = snakemake.input.library
    output_csv = snakemake.output.csv
    output_smi = snakemake.output.smi

    require_file(scores_path, "tabla de scores")
    require_file(library_path, "SMILES enumerados")
    ensure_parent_dir(output_csv)
    ensure_parent_dir(output_smi)

    top_percent = clamp_top_percent(getattr(snakemake.params, "top_percent", 5))

    scores_df = pd.read_csv(scores_path)
    if scores_df.empty:
        log_warn("⚠️ Tabla de scores vacía. No se pueden seleccionar candidatos.")
        pd.DataFrame(columns=["ligand_id", "affinity", "smiles", "original_id"]).to_csv(
            output_csv, index=False
        )
        with open(output_smi, "w", encoding="utf-8") as handle:
            handle.write("")
        confirm_file(output_csv, "candidatos top")
        confirm_file(output_smi, "candidatos top (smi)")
        return

    if not {"ligand_id", "affinity"}.issubset(scores_df.columns):
        log_warn("⚠️ El CSV de scores no tiene columnas requeridas: ligand_id, affinity")
        scores_df = scores_df.rename(
            columns={scores_df.columns[0]: "ligand_id", scores_df.columns[1]: "affinity"}
        )

    scores_df = scores_df.sort_values("affinity", ascending=True)
    total = len(scores_df)
    top_count = max(1, math.ceil(total * top_percent / 100))

    log_info(
        "✅ Seleccionando top "
        f"{top_percent:.1f}% ({top_count}/{total}) por afinidad."
    )

    top_df = scores_df.head(top_count).copy()

    smiles_map = load_library_smi(library_path)

    smiles_list = []
    original_ids = []
    missing = 0
    for ligand_id in top_df["ligand_id"]:
        match = smiles_map.get(ligand_id)
        if match is None:
            missing += 1
            smiles_list.append("")
            original_ids.append("")
        else:
            smiles, original_id = match
            smiles_list.append(smiles)
            original_ids.append(original_id)

    if missing:
        log_warn(f"⚠️ {missing} candidatos sin SMILES asociado en la librería.")

    top_df["smiles"] = smiles_list
    top_df["original_id"] = original_ids
    top_df.to_csv(output_csv, index=False)

    with open(output_smi, "w", encoding="utf-8") as handle:
        for ligand_id, smiles in zip(top_df["ligand_id"], top_df["smiles"]):
            if not smiles:
                continue
            handle.write(f"{smiles}\t{ligand_id}\n")

    confirm_file(output_csv, "candidatos top")
    confirm_file(output_smi, "candidatos top (smi)")
    log_info(f"📦 Output listo para Vina/Gnina: {output_smi}")


if __name__ == "__main__":
    main()
