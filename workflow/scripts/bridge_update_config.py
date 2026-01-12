# --- MOCK PARA DESARROLLO (Pylance no se quejará) ---
if "snakemake" not in globals():
    # Esto simula la variable snakemake para que el editor no marque error
    # y tú sepas qué estructura tiene. Solo sirve mientras editas.
    from types import SimpleNamespace
    snakemake = SimpleNamespace(
        input=SimpleNamespace(pdb="", fasta="", xml="", original_fasta="", msa=""),
        output=SimpleNamespace(pdb="", fasta="", xml="", fasta_msa="", pdb_conserved=""),
        params=SimpleNamespace(pdb_id="7KGY", chain="B", n_hits=10, e_val=0.001),
        wildcards=SimpleNamespace()
    )
# ----------------------------------------------------
import pandas as pd
import yaml
import numpy as np
from Bio.PDB import PDBParser
import sys

# Inputs desde Snakemake
ranking_csv = snakemake.input.csv
pdb_file = snakemake.input.pdb
config_file = "config/config.yaml"

print("🌉 PUENTE DE DATOS: Actualizando configuración con el mejor blanco descubierto...")

# 1. Leer el Ganador
try:
    df = pd.read_csv(ranking_csv)
    # El mejor es la primera fila (ya está ordenado por Ghost_Score)
    winner = df.iloc[0]
    best_res_id = int(winner['ID'])
    best_score = winner['Ghost_Score']
    
    print(f"   🏆 Ganador Indiscutible: Residuo {winner['Residue']}{best_res_id}")
    print(f"   ⭐ Ghost Score: {best_score}")

except Exception as e:
    print(f"❌ Error leyendo el ranking: {e}")
    sys.exit(1)

# 2. Calcular Centro Geométrico (X, Y, Z) del Ganador
parser = PDBParser(QUIET=True)
structure = parser.get_structure("Target", pdb_file)
chain_id = snakemake.params.chain

target_atoms = []
for model in structure:
    for chain in model:
        if chain.id == chain_id:
            for residue in chain:
                if residue.id[1] == best_res_id:
                    target_atoms = [atom.get_coord() for atom in residue]
                    break

if not target_atoms:
    print(f"❌ Error: El residuo {best_res_id} no tiene átomos en el PDB.")
    sys.exit(1)

center = np.mean(target_atoms, axis=0)
cx, cy, cz = float(center[0]), float(center[1]), float(center[2])

print(f"   📍 Coordenadas calculadas: [{cx:.2f}, {cy:.2f}, {cz:.2f}]")

# 3. Reescribir config.yaml
with open(config_file, 'r') as f:
    config = yaml.safe_load(f)

# Actualizamos valores
config['structure']['target_residue'] = best_res_id
config['docking']['center_x'] = cx
config['docking']['center_y'] = cy
config['docking']['center_z'] = cz

# Escribimos de vuelta preservando el orden (best effort)
with open(config_file, 'w') as f:
    yaml.dump(config, f, sort_keys=False)

print("✅ Configuración actualizada automáticamente. El Módulo 2 apuntará al lugar correcto.")