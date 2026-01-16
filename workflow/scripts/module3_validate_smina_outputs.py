# --- MOCK PARA DESARROLLO (Pylance no se quejará) ---
if "snakemake" not in globals():
    from types import SimpleNamespace

    snakemake = SimpleNamespace(
        input=SimpleNamespace(docked=[]),
        output=SimpleNamespace(flag="results/module3/smina_outputs_validated.flag"),
    )

from pathlib import Path

from logging_utils import confirm_file, ensure_parent_dir, log_error, log_info


def main():
    docked_inputs = list(getattr(snakemake.input, "docked", []))
    docked_paths = [Path(path) for path in docked_inputs if path]
    existing = [path for path in docked_paths if path.is_file()]

    log_info(
        "🧪 Validando poses de Smina en results/module3/docking: "
        f"{len(existing)} de {len(docked_paths)} archivos encontrados."
    )

    if not existing:
        log_error(
            "No se encontraron poses de Smina (*.pdbqt) en results/module3/docking. "
            "Deteniendo el flujo para evitar tablas vacías."
        )
        raise SystemExit(1)

    output_flag = getattr(snakemake.output, "flag", None) or snakemake.output[0]
    ensure_parent_dir(output_flag)
    Path(output_flag).write_text("validated\n", encoding="utf-8")
    confirm_file(output_flag, "validación de poses de Smina")
    log_info(f"✅ Validación completada: {output_flag}")


if __name__ == "__main__":
    main()
