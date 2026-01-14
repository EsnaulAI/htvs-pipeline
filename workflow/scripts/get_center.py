# workflow/scripts/get_center.py

# --- MOCK PARA DESARROLLO ---
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
        input=SimpleNamespace(pdb=config["evolution"]["conservation_pdb"]),
        output=SimpleNamespace(coords=config["docking"]["gridbox_center"]),
        params=SimpleNamespace(
            chain=config["structure"]["chain_id"],
            target_res=config["structure"]["target_residue"],
        )
    )
# ----------------------------

from pathlib import Path
import numpy as np
from Bio.PDB import PDBParser
import sys
import yaml

def load_config():
    config_path = Path(__file__).resolve().parents[2] / "config" / "config.yaml"
    with open(config_path, "r") as f:
        return yaml.safe_load(f)

pdb_file = snakemake.input.pdb
chain_id = snakemake.params.chain
target_res_num = int(snakemake.params.target_res)
output_file = snakemake.output.coords

print(f"📍 Calculando centroide geométrico para {chain_id}:{target_res_num}...")

parser = PDBParser(QUIET=True)
structure = parser.get_structure('target', pdb_file)

target_residue = None
found = False

# Buscar el residuo
for model in structure:
    for chain in model:
        if chain.id == chain_id:
            for residue in chain:
                if residue.id[1] == target_res_num:
                    target_residue = residue
                    found = True
                    break
        if found: break
    if found: break

if not target_residue:
    print(f"❌ Error: Residuo {target_res_num} no encontrado en cadena {chain_id}.")
    sys.exit(1)

# Calcular centroide
coords = []
for atom in target_residue:
    coords.append(atom.get_coord())

center = np.mean(coords, axis=0)
x, y, z = center[0], center[1], center[2]

print(f"✅ Centroide calculado: [{x:.3f}, {y:.3f}, {z:.3f}]")

# Guardar
with open(output_file, "w") as f:
    f.write(f"center_x = {x:.3f}\n")
    f.write(f"center_y = {y:.3f}\n")
    f.write(f"center_z = {z:.3f}\n")
