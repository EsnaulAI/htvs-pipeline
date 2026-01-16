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
import networkx as nx
from prody import parsePDB, ANM, calcSqFlucts, confProDy
from Bio.PDB import PDBParser, SASA
import numpy as np
import re
import sys
import os
import glob
from logging_utils import (
    confirm_file,
    ensure_parent_dir,
    log_info,
    log_error,
    require_file,
)

# Inputs
pdb_file = snakemake.input.pdb
pocket_file = snakemake.input.pockets
chain_id = snakemake.params.chain
output_csv = snakemake.output.csv

confProDy(verbosity='none')

def get_pocket_residues(file_path):
    """Parsea el output de fpocket o un PDB de cavidades para extraer residuos en bolsillos."""
    pocket_res = set()
    if not file_path or not os.path.exists(file_path):
        return pocket_res

    def parse_pocket_pdb(pdb_path):
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("Pocket", pdb_path)
        for model in structure:
            for chain in model:
                for res in chain:
                    if res.id[0] == " ":
                        pocket_res.add(res.id[1])

    if file_path.lower().endswith(".pdb"):
        parse_pocket_pdb(file_path)
        return pocket_res

    try:
        with open(file_path, "r") as f:
            for line in f:
                if re.search(r"\bResidues?\b", line, re.IGNORECASE):
                    pocket_res.update(int(res_id) for res_id in re.findall(r"\b\d+\b", line))
    except OSError:
        return pocket_res

    if not pocket_res:
        base_dir = os.path.dirname(file_path)
        candidates = sorted(
            glob.glob(os.path.join(base_dir, "*pocket*.pdb")) +
            glob.glob(os.path.join(base_dir, "*pocket*_atm.pdb")) +
            glob.glob(os.path.join(base_dir, "*pocket*_vert.pdb")) +
            glob.glob(os.path.join(base_dir, "pocket*.pdb"))
        )
        for candidate in candidates:
            parse_pocket_pdb(candidate)
    return pocket_res

require_file(pdb_file, "PDB conservado")
require_file(pocket_file, "archivo de pockets")
ensure_parent_dir(output_csv)

log_info(f"🕵️ Explorando candidatos en cadena {chain_id}...")

# Residuos en cavidades reales (fpocket)
pocket_residues = get_pocket_residues(pocket_file)
print(f"   🧩 Residuos en pockets detectados: {len(pocket_residues)}")

# 1. Cargar Estructura
parser = PDBParser(QUIET=True)
structure = parser.get_structure("Target", pdb_file)
chain = structure[0][chain_id]

# 2. Calcular SASA (Accesibilidad) - CRÍTICO PARA NO ELEGIR SITIOS ENTERRADOS
sr = SASA.ShrakeRupley()
sr.compute(structure, level="R")

# 3. Network Analysis
G = nx.Graph()
residues = [r for r in chain if r.id[0] == " "]
for i, r1 in enumerate(residues):
    for j, r2 in enumerate(residues):
        if i >= j: continue
        try:
            if (r1['CA'] - r2['CA']) < 8.0:
                G.add_edge(r1.id[1], r2.id[1])
        except: continue
centrality = nx.betweenness_centrality(G)

# 4. Dinámica (NMA - Rigidez)
prody_struct = parsePDB(pdb_file)
calphas = prody_struct.select(f'chain {chain_id} and calpha')
if calphas is None or calphas.numAtoms() == 0:
    log_error(
        f"No se encontraron C-alfa para cadena {chain_id}; no se puede construir ANM."
    )
    sys.exit(1)
anm = ANM('Analysis')
anm.buildHessian(calphas)
anm.calcModes(n_modes=20)
mobility_map = dict(zip(calphas.getResnums(), calcSqFlucts(anm)))

# 5. INTEGRACIÓN (Ghost-Score v2.0)
data = []

log_info("📊 Evaluando residuos...")

sasa_values = [res.sasa if hasattr(res, 'sasa') else 0 for res in residues]
sasa_min = min(sasa_values) if sasa_values else 0
sasa_max = max(sasa_values) if sasa_values else 0
sasa_range = max(sasa_max - sasa_min, 1e-6)

for res in residues:
    res_id = res.id[1]
    
    # Métricas Base
    cons_score = res['CA'].get_bfactor() if 'CA' in res else 0 # Viene de map_conservation
    cent_score = centrality.get(res_id, 0)
    
    mob = mobility_map.get(res_id, 999)
    rigidity = 1 / mob if mob > 0 else 0
    
    sasa = res.sasa if hasattr(res, 'sasa') else 0
    
    # --- FILTRO DE DRUGGABILITY ---
    sasa_norm = (sasa - sasa_min) / sasa_range
    pocket_member = 1.0 if res_id in pocket_residues else 0.0
    # Fórmula nueva: pocket_term = (Wpocket * pocket_member) + (Wsasa * SASA_norm)
    pocket_weight = 20.0
    sasa_weight = 10.0
    pocket_term = (pocket_weight * pocket_member) + (sasa_weight * sasa_norm)
        
    # Fórmula Maestra
    final_score = (cons_score * 3.0) + (cent_score * 50.0) + (rigidity * 1.5) + pocket_term
    
    data.append({
        "Residue": res.get_resname(),
        "ID": res_id,
        "Ghost_Score": round(final_score, 3),
        "Conservation": round(cons_score, 2),
        "SASA": round(sasa, 1),
        "SASA_Norm": round(sasa_norm, 3),
        "Pocket_Member": int(pocket_member),
        "Pocket_Term": round(pocket_term, 3),
        "Rigidity": round(rigidity, 1)
    })

# Ranking
df = pd.DataFrame(data).sort_values(by="Ghost_Score", ascending=False)
df.to_csv(output_csv, index=False)
confirm_file(output_csv, "ranking de candidatos")

winner = df.iloc[0]
print(f"🏆 GANADOR CIENTÍFICO: {winner['Residue']}{winner['ID']}")
print(f"   Score: {winner['Ghost_Score']} (SASA: {winner['SASA']})")
log_info(f"🏆 GANADOR CIENTÍFICO: {winner['Residue']}{winner['ID']}")
log_info(f"Score: {winner['Ghost_Score']} (SASA: {winner['SASA']})")
