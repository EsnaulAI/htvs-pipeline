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
        input=SimpleNamespace(gnina_scores="results/module3/gnina_scores.csv"),
        output=SimpleNamespace(consensus="results/module3/consensus_final.csv"),
        wildcards=SimpleNamespace(),
    )

import pandas as pd

from logging_utils import (
    confirm_file,
    ensure_parent_dir,
    log_info,
    log_warn,
    require_file,
)


def main():
    input_scores = snakemake.input.gnina_scores
    output_csv = snakemake.output.consensus

    require_file(input_scores, "scores Gnina")
    ensure_parent_dir(output_csv)

    scores_df = pd.read_csv(input_scores)
    if scores_df.empty:
        log_warn("⚠️ Tabla Gnina vacía. Consenso sin datos.")
        pd.DataFrame(
            columns=[
                "ligand_id",
                "vina_score",
                "cnn_score",
                "vina_rank",
                "cnn_rank",
                "consensus_rank",
            ]
        ).to_csv(output_csv, index=False)
        confirm_file(output_csv, "consenso final")
        return

    if not {"ligand_id", "vina_score", "cnn_score"}.issubset(scores_df.columns):
        log_warn("⚠️ Columnas faltantes en gnina_scores.csv; usando primeras columnas.")
        scores_df = scores_df.rename(
            columns={
                scores_df.columns[0]: "ligand_id",
                scores_df.columns[1]: "vina_score",
                scores_df.columns[2]: "cnn_score",
            }
        )

    scores_df["vina_rank"] = scores_df["vina_score"].rank(
        ascending=True, method="min"
    )
    scores_df["cnn_rank"] = scores_df["cnn_score"].rank(
        ascending=False, method="min"
    )
    scores_df["consensus_rank"] = scores_df["vina_rank"] + scores_df["cnn_rank"]

    scores_df = scores_df.sort_values("consensus_rank", ascending=True)
    scores_df.to_csv(output_csv, index=False)
    confirm_file(output_csv, "consenso final")
    log_info(f"✅ Consenso final listo: {output_csv}")


if __name__ == "__main__":
    main()
