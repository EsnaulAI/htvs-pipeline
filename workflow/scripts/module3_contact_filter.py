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
            receptor=config["evolution"]["conservation_pdb"],
        ),
        output=SimpleNamespace(final_hits="results/module3/final_hits.csv"),
        params=SimpleNamespace(
            chain=config["structure"]["chain_id"],
            target_res=config["structure"]["target_residue"],
            pocket_radius=6.0,
            contact_cutoff=4.0,
            min_contacts=2,
            require_target_contact=True,
        ),
        wildcards=SimpleNamespace(),
    )

import os
from typing import Dict, Iterable, List, Tuple

import numpy as np
import pandas as pd
from Bio.PDB import NeighborSearch, PDBParser

from logging_utils import (
    confirm_file,
    ensure_parent_dir,
    log_info,
    log_warn,
    require_file,
)


def parse_pdbqt_atoms(pose_path: str) -> np.ndarray:
    coords: List[Tuple[float, float, float]] = []
    with open(pose_path, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                except ValueError:
                    parts = line.split()
                    if len(parts) < 9:
                        continue
                    x, y, z = map(float, parts[6:9])
                coords.append((x, y, z))
    if not coords:
        return np.empty((0, 3))
    return np.array(coords, dtype=float)


def build_key_residues(
    receptor_path: str,
    chain_id: str,
    target_res: int,
    pocket_radius: float,
) -> Dict[int, np.ndarray]:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("receptor", receptor_path)
    chain = structure[0][chain_id]

    residues = [res for res in chain.get_residues() if res.id[0] == " "]
    target_residue = next((res for res in residues if res.id[1] == target_res), None)
    if target_residue is None:
        raise ValueError(f"Residuo objetivo {target_res} no encontrado en {chain_id}.")

    center_atom = target_residue["CA"] if "CA" in target_residue else target_residue.child_list[0]
    atom_list = [atom for atom in chain.get_atoms()]
    ns = NeighborSearch(atom_list)
    neighbors = ns.search(center_atom.get_coord(), pocket_radius, level="R")
    key_residues = {res.id[1]: res for res in neighbors if res.id[0] == " "}
    key_residues[target_res] = target_residue

    residue_coords: Dict[int, np.ndarray] = {}
    for res_id, residue in key_residues.items():
        atoms = [atom.get_coord() for atom in residue.get_atoms()]
        if atoms:
            residue_coords[res_id] = np.array(atoms, dtype=float)
    return residue_coords


def residue_contacts(
    ligand_coords: np.ndarray,
    residue_coords: np.ndarray,
    cutoff: float,
) -> bool:
    if ligand_coords.size == 0 or residue_coords.size == 0:
        return False
    diff = ligand_coords[:, None, :] - residue_coords[None, :, :]
    dist2 = np.sum(diff * diff, axis=2)
    return np.any(dist2 <= cutoff * cutoff)


def collect_pose_rows(scores_path: str) -> pd.DataFrame:
    scores_df = pd.read_csv(scores_path)
    if scores_df.empty:
        return pd.DataFrame(columns=["ligand_id", "score", "pose_path"])
    if "pose_path" not in scores_df.columns:
        if {"ligand_id", "score"}.issubset(scores_df.columns):
            scores_df = scores_df.assign(
                pose_path=scores_df["ligand_id"].apply(
                    lambda x: os.path.join("results/module3/docking", f"{x}_out.pdbqt")
                )
            )
        else:
            scores_df = scores_df.rename(
                columns={
                    scores_df.columns[0]: "ligand_id",
                    scores_df.columns[1]: "score",
                }
            )
            scores_df = scores_df.assign(
                pose_path=scores_df["ligand_id"].apply(
                    lambda x: os.path.join("results/module3/docking", f"{x}_out.pdbqt")
                )
            )
    return scores_df


def format_residue_list(res_ids: Iterable[int]) -> str:
    return ",".join(str(res_id) for res_id in sorted(res_ids))


def main():
    scores_path = snakemake.input.scores
    receptor_path = snakemake.input.receptor
    output_csv = snakemake.output.final_hits

    chain_id = snakemake.params.chain
    target_res = int(snakemake.params.target_res)
    pocket_radius = float(getattr(snakemake.params, "pocket_radius", 6.0))
    contact_cutoff = float(getattr(snakemake.params, "contact_cutoff", 4.0))
    min_contacts = int(getattr(snakemake.params, "min_contacts", 2))
    require_target_contact = bool(getattr(snakemake.params, "require_target_contact", True))

    require_file(scores_path, "scores de docking")
    require_file(receptor_path, "receptor PDB")
    ensure_parent_dir(output_csv)

    log_info("🔍 Construyendo residuos clave para filtrado de contactos...")
    try:
        residue_coords = build_key_residues(
            receptor_path=receptor_path,
            chain_id=chain_id,
            target_res=target_res,
            pocket_radius=pocket_radius,
        )
    except ValueError as exc:
        log_warn(f"⚠️ {exc}")
        residue_coords = {}

    if not residue_coords:
        log_warn("⚠️ No se encontraron residuos clave; se marcarán sin contactos.")

    key_residue_ids = sorted(residue_coords.keys())
    log_info(
        "📌 Residuos clave: "
        f"{format_residue_list(key_residue_ids)} (total={len(key_residue_ids)})"
    )

    scores_df = collect_pose_rows(scores_path)
    if scores_df.empty:
        log_warn("⚠️ Tabla de docking vacía; no hay poses para filtrar.")
        pd.DataFrame(
            columns=[
                "ligand_id",
                "score",
                "pose_path",
                "contacts",
                "target_contact",
                "contact_residues",
                "pass_contact_filter",
            ]
        ).to_csv(output_csv, index=False)
        confirm_file(output_csv, "final hits")
        return

    records = []
    for _, row in scores_df.iterrows():
        ligand_id = str(row.get("ligand_id", "")).strip()
        pose_path = row.get("pose_path", "")
        score = row.get("score", row.get("affinity", None))

        if not pose_path or not os.path.isfile(pose_path):
            log_warn(f"⚠️ Pose inexistente para {ligand_id}: {pose_path}")
            contacts = 0
            contacted_residues: List[int] = []
        else:
            ligand_coords = parse_pdbqt_atoms(pose_path)
            contacted_residues = []
            for res_id, coords in residue_coords.items():
                if residue_contacts(ligand_coords, coords, contact_cutoff):
                    contacted_residues.append(res_id)
            contacts = len(contacted_residues)

        target_contact = target_res in contacted_residues
        passes_contacts = contacts >= min_contacts
        if require_target_contact:
            passes_contacts = passes_contacts and target_contact

        records.append(
            {
                "ligand_id": ligand_id,
                "score": score,
                "pose_path": pose_path,
                "contacts": contacts,
                "target_contact": target_contact,
                "contact_residues": format_residue_list(contacted_residues),
                "pass_contact_filter": passes_contacts,
            }
        )

    df = pd.DataFrame(records)
    final_hits = df[df["pass_contact_filter"]].copy()
    if final_hits.empty:
        log_warn("⚠️ Ningún ligando cumplió el filtro de contactos.")
    final_hits.to_csv(output_csv, index=False)

    confirm_file(output_csv, "final hits")
    log_info(f"✅ Filtro de contactos completado: {len(final_hits)}/{len(df)} hits.")


if __name__ == "__main__":
    main()
