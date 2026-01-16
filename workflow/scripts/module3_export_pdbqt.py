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

INVALID_LIGANDS_LOG = os.path.join("results", "module3", "invalid_ligands.txt")

def has_atoms(pdbqt_path: str) -> bool:
    if not os.path.isfile(pdbqt_path):
        return False
    with open(pdbqt_path, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if line.startswith(("ATOM", "HETATM")):
                return True
    return False


def has_pdbqt_atoms(pdbqt_path):
    if not os.path.isfile(pdbqt_path):
        return False
    with open(pdbqt_path, "r", encoding="utf-8", errors="ignore") as handle:
        return any(
            line.startswith("ATOM") or line.startswith("HETATM")
            for line in handle
        )


def log_invalid_ligands(invalid_ids, log_path):
    os.makedirs(os.path.dirname(log_path), exist_ok=True)
    existing = set()
    if os.path.isfile(log_path):
        with open(log_path, "r", encoding="utf-8") as handle:
            existing = {line.strip() for line in handle if line.strip()}
    updated = list(existing.union(invalid_ids))
    with open(log_path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(sorted(updated)))


def filter_ligands_list(ligands_list_path, invalid_ids):
    if not os.path.isfile(ligands_list_path):
        return False
    with open(ligands_list_path, "r", encoding="utf-8") as handle:
        ligands = [line.strip() for line in handle if line.strip()]
    remaining = [ligand for ligand in ligands if ligand not in invalid_ids]
    with open(ligands_list_path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(remaining))
    return True


def convert_sdf_to_pdbqt(sdf_path, out_dir, invalid_log_path):
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
    if not has_atoms(out_path):
        log_warn(f"PDBQT vacío o sin átomos: {out_path}")
        try:
            os.remove(out_path)
        except FileNotFoundError:
            pass
        return False
    return True
        return False, None
    if not os.path.isfile(out_path) or os.path.getsize(out_path) == 0:
        log_warn(f"PDBQT vacío o inexistente: {out_path}")
        if os.path.isfile(out_path):
            os.remove(out_path)
        log_invalid_ligands({safe_name}, invalid_log_path)
        return False, safe_name
    if not has_pdbqt_atoms(out_path):
        log_warn(f"PDBQT sin átomos para {safe_name}; eliminado.")
        os.remove(out_path)
        log_invalid_ligands({safe_name}, invalid_log_path)
        return False, safe_name
    return True, None


def pdbqt_has_atoms(pdbqt_path):
    with open(pdbqt_path, "r", encoding="utf-8") as pdbqt_file:
        for line in pdbqt_file:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                return True
    return False


def main():
    sdf_dir = snakemake.input.sdf_dir
    out_dir = snakemake.output.pdbqt_dir

    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(os.path.dirname(INVALID_LIGANDS_LOG), exist_ok=True)
    invalid_log_path = os.path.join("results", "module3", "invalid_ligands.txt")

    sdf_files = sorted(glob.glob(os.path.join(sdf_dir, "*.sdf")))
    if not sdf_files:
        log_warn(f"No se encontraron ligandos SDF en {sdf_dir}.")

    converted = 0
    valid = 0
    invalid = 0
    for sdf_path in sdf_files:
        converted_ok = convert_sdf_to_pdbqt(sdf_path, out_dir)
        if converted_ok:
            converted += 1
            name = os.path.splitext(os.path.basename(sdf_path))[0]
            safe_name = "".join([c for c in name if c.isalnum() or c in ("_", "-")])
            if not safe_name:
                safe_name = name
            out_path = os.path.join(out_dir, f"{safe_name}.pdbqt")
            if not pdbqt_has_atoms(out_path):
                os.remove(out_path)
                with open(INVALID_LIGANDS_LOG, "a", encoding="utf-8") as log_file:
                    log_file.write(f"{safe_name}\n")
                log_warn(f"PDBQT sin átomos, eliminado: {out_path}")
                invalid += 1
            else:
                valid += 1
    invalid_ids = set()
    for sdf_path in sdf_files:
        success, invalid_id = convert_sdf_to_pdbqt(
            sdf_path, out_dir, invalid_log_path
        )
        if success:
            converted += 1
        elif invalid_id:
            invalid_ids.add(invalid_id)

    if invalid_ids:
        if filter_ligands_list("results/module2/ligands_list.txt", invalid_ids):
            log_info(
                "🧹 Ligandos inválidos eliminados del listado para docking: "
                f"{len(invalid_ids)}"
            )
        log_info(
            "📄 Ligandos inválidos registrados en "
            f"{invalid_log_path}: {len(invalid_ids)}"
        )

    confirm_path(out_dir, "directorio de salida PDBQT")
    log_info(f"✅ Exportación completada: {converted} ligandos en {out_dir}")
    log_info(f"Resumen: válidos={valid}, inválidos={invalid}")


if __name__ == "__main__":
    main()
