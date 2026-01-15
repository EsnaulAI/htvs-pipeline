# --- MOCK PARA DESARROLLO (Pylance no se quejará) ---
if "snakemake" not in globals():
    from types import SimpleNamespace
    snakemake = SimpleNamespace(
        input=SimpleNamespace(sdf_dir=""),
        output=SimpleNamespace(pdbqt_dir=""),
        params=SimpleNamespace(),
        wildcards=SimpleNamespace()
    )
# ----------------------------------------------------
# workflow/scripts/module3_export_pdbqt.py
import glob
import os
import subprocess

from logging_utils import (
    confirm_path,
    log_info,
    log_warn,
)


def convert_sdf_to_pdbqt(sdf_path, out_dir):
    name = os.path.splitext(os.path.basename(sdf_path))[0]
    safe_name = "".join([c for c in name if c.isalnum() or c in ("_", "-")])
    if not safe_name:
        safe_name = name
    out_path = os.path.join(out_dir, f"{safe_name}.pdbqt")
    command = [
        "obabel",
        sdf_path,
        "-O",
        out_path,
        "--partialcharge",
        "gasteiger",
    ]
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        log_warn(
            "Falló obabel para {name}: {stderr}".format(
                name=safe_name,
                stderr=result.stderr.strip() or "sin salida de error",
            )
        )
        return False
    if not os.path.isfile(out_path) or os.path.getsize(out_path) == 0:
        log_warn(f"PDBQT vacío o inexistente: {out_path}")
        return False
    return True


def main():
    sdf_dir = snakemake.input.sdf_dir
    out_dir = snakemake.output.pdbqt_dir

    os.makedirs(out_dir, exist_ok=True)

    sdf_files = sorted(glob.glob(os.path.join(sdf_dir, "*.sdf")))
    if not sdf_files:
        log_warn(f"No se encontraron ligandos SDF en {sdf_dir}.")

    converted = 0
    for sdf_path in sdf_files:
        if convert_sdf_to_pdbqt(sdf_path, out_dir):
            converted += 1

    confirm_path(out_dir, "directorio de salida PDBQT")
    log_info(f"✅ Exportación completada: {converted} ligandos en {out_dir}")


if __name__ == "__main__":
    main()
