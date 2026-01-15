from __future__ import annotations

import csv
import subprocess
from pathlib import Path

receptor = Path(snakemake.input.receptor)
top_candidates = Path(snakemake.input.top_candidates)
ligands_dir = Path(snakemake.params.ligands_dir)
vina_dir = Path(snakemake.output.vina_dir)
vina_logs = Path(snakemake.output.vina_logs)

vina_dir.mkdir(parents=True, exist_ok=True)
vina_logs.mkdir(parents=True, exist_ok=True)

center_x = str(snakemake.params.center_x)
center_y = str(snakemake.params.center_y)
center_z = str(snakemake.params.center_z)
size_x = str(snakemake.params.size_x)
size_y = str(snakemake.params.size_y)
size_z = str(snakemake.params.size_z)

with top_candidates.open(newline="", encoding="utf-8") as handle:
    reader = csv.DictReader(handle)
    for row in reader:
        ligand_id = row.get("ligand_id") or row.get("ligand") or ""
        ligand_id = ligand_id.strip()
        if not ligand_id:
            continue
        ligand_path = ligands_dir / f"{ligand_id}.pdbqt"
        if not ligand_path.exists():
            continue
        out_path = vina_dir / f"{ligand_id}_vina_out.pdbqt"
        log_path = vina_logs / f"{ligand_id}.log"
        cmd = [
            "vina",
            "--receptor",
            str(receptor),
            "--ligand",
            str(ligand_path),
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
            str(out_path),
            "--log",
            str(log_path),
        ]
        subprocess.run(cmd, check=True)
