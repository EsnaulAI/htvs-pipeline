# --- MOCK PARA DESARROLLO (Pylance no se quejará) ---
if "snakemake" not in globals():
    from types import SimpleNamespace
    snakemake = SimpleNamespace(
        input=SimpleNamespace(prepared_done="", prepared_dir=""),
        output=SimpleNamespace(pdbqt_dir=""),
        params=SimpleNamespace(),
        wildcards=SimpleNamespace()
    )
# ----------------------------------------------------
# workflow/scripts/module2_export_pdbqt.py
import glob
import os
from openbabel import pybel
from logging_utils import (
    confirm_path,
    log_info,
    log_warn,
    require_file,
)


def export_pdbqt(sdf_path, out_dir):
    for mol in pybel.readfile("sdf", sdf_path):
        name = mol.title or os.path.splitext(os.path.basename(sdf_path))[0]
        safe_name = "".join([c for c in name if c.isalnum() or c in ("_", "-")])
        if not safe_name:
            safe_name = os.path.splitext(os.path.basename(sdf_path))[0]
        out_path = os.path.join(out_dir, f"{safe_name}.pdbqt")
        mol.addh()
        try:
            mol.calccharges("gasteiger")
        except Exception:
            log_warn(f"No se pudieron calcular cargas para {safe_name}; exportando sin cargas.")
        mol.write("pdbqt", out_path, overwrite=True)


def main():
    prepared_done = snakemake.input.prepared_done
    sdf_dir = snakemake.input.prepared_dir
    out_dir = snakemake.output.pdbqt_dir

    require_file(prepared_done, "flag de preparación 3D")

    os.makedirs(out_dir, exist_ok=True)

    sdf_files = sorted(glob.glob(os.path.join(sdf_dir, "*.sdf")))
    if not sdf_files:
        log_warn(f"No se encontraron ligandos SDF en {sdf_dir}.")

    count = 0
    for sdf_path in sdf_files:
        export_pdbqt(sdf_path, out_dir)
        count += 1

    confirm_path(out_dir, "directorio de salida PDBQT")
    log_info(f"✅ Exportación completada: {count} ligandos en {out_dir}")


if __name__ == "__main__":
    main()
