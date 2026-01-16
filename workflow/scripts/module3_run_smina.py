from __future__ import annotations

from pathlib import Path
import subprocess

from logging_utils import (
    ensure_parent_dir,
    log_info,
    log_warn,
    require_file,
)


def ligand_has_atoms(ligand_path: Path) -> bool:
    if not ligand_path.is_file():
        return False
    with ligand_path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if line.startswith(("ATOM", "HETATM")):
                return True
    return False


def main() -> None:
    receptor = Path(snakemake.input.receptor)
    ligand = Path(snakemake.input.ligand)
    docked = Path(snakemake.output.docked)
    log_path = Path(snakemake.output.log)

    require_file(receptor, "receptor preparado")
    require_file(ligand, "ligando PDBQT")
    ensure_parent_dir(docked)
    ensure_parent_dir(log_path)

    if not ligand_has_atoms(ligand):
        log_warn(f"⚠️ Ligando sin átomos: {ligand}. Se omite docking.")
        log_path.write_text("SKIPPED: ligand has no atoms\n", encoding="utf-8")
        docked.write_text("", encoding="utf-8")
        return

import subprocess
from pathlib import Path

from logging_utils import confirm_file, ensure_parent_dir, log_info, log_warn, require_file


receptor = Path(snakemake.input.receptor)
ligand = Path(snakemake.input.ligand)
docked = Path(snakemake.output.docked)
log_path = Path(snakemake.output.log)

center_x = str(snakemake.params.cx)
center_y = str(snakemake.params.cy)
center_z = str(snakemake.params.cz)
size_x = str(snakemake.params.sx)
size_y = str(snakemake.params.sy)
size_z = str(snakemake.params.sz)

require_file(str(receptor), label="receptor PDBQT")
require_file(str(ligand), label="ligando PDBQT")

ensure_parent_dir(str(docked))
ensure_parent_dir(str(log_path))

has_atoms = False
with ligand.open("r", encoding="utf-8") as handle:
    for line in handle:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            has_atoms = True
            break

if not has_atoms:
    log_warn(f"Ligando sin átomos en {ligand}. Se omite Smina.")
    docked.write_text("", encoding="utf-8")
    log_path.write_text("SKIPPED: ligando sin átomos\n", encoding="utf-8")
else:
    cmd = [
        "smina",
        "--receptor",
        str(receptor),
        "--ligand",
        str(ligand),
        "--center_x",
        str(snakemake.params.cx),
        "--center_y",
        str(snakemake.params.cy),
        "--center_z",
        str(snakemake.params.cz),
        "--size_x",
        str(snakemake.params.sx),
        "--size_y",
        str(snakemake.params.sy),
        "--size_z",
        str(snakemake.params.sz),
        center_x,
        "--center_y",
        center_y,
        "--center_z",
        center_z,
        "--size_x",
        size_x,
        "--size_y",
        size_y,
        "--size_z",
        size_z,
        "--out",
        str(docked),
        "--log",
        str(log_path),
        "--cpu",
        "1",
    ]
    log_info(f"🧪 Ejecutando Smina para {ligand.name}...")
    subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()
    log_info("Ejecutando Smina.")
    subprocess.run(cmd, check=True)
    confirm_file(str(docked), label="pose de docking")
    confirm_file(str(log_path), label="log de Smina")
