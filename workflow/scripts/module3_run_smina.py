from __future__ import annotations

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
    log_info("Ejecutando Smina.")
    subprocess.run(cmd, check=True)
    confirm_file(str(docked), label="pose de docking")
    confirm_file(str(log_path), label="log de Smina")
