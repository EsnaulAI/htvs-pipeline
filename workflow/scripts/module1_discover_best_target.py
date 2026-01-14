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
from logging_utils import (
    confirm_file,
    ensure_parent_dir,
    log_info,
    require_file,
)

# Inputs
pdb_file = snakemake.input.pdb
pocket_file = snakemake.input.pockets
chain_id = snakemake.params.chain
output_csv = snakemake.output.csv

confProDy(verbosity='none')

def get_pocket_residues(file_path):
    """Parsea el archivo info de fpocket para sacar residuos en cavidades."""
    pocket_res = set()
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()
            # Buscamos patrones simples de índices de residuos si existen en el reporte
            # Nota: fpocket info es complejo. Una heurística simple es confiar 
            # en la combinación de SASA y el archivo PDB de pockets si se usara.
            # Aquí usamos SASA como proxy principal y fpocket como bono si logramos parsearlo.
            pass 
    except:
        pass
    return pocket_res

require_file(pdb_file, "PDB conservado")
require_file(pocket_file, "archivo de pockets")
ensure_parent_dir(output_csv)

log_info(f"🕵️ Explorando candidatos en cadena {chain_id}...")

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
anm = ANM('Analysis')
anm.buildHessian(calphas)
anm.calcModes(n_modes=20)
mobility_map = dict(zip(calphas.getResnums(), calcSqFlucts(anm)))

# 5. INTEGRACIÓN (Ghost-Score v2.0)
data = []

log_info("📊 Evaluando residuos...")

for res in residues:
    res_id = res.id[1]
    
    # Métricas Base
    cons_score = res['CA'].get_bfactor() if 'CA' in res else 0 # Viene de map_conservation
    cent_score = centrality.get(res_id, 0)
    
    mob = mobility_map.get(res_id, 999)
    rigidity = 1 / mob if mob > 0 else 0
    
    sasa = res.sasa if hasattr(res, 'sasa') else 0
    
    # --- FILTRO DE DRUGGABILITY ---
    pocket_bonus = 0.0
    
    # Rango ideal de SASA para un bolsillo: 10 - 80 Angstroms^2
    if 10.0 <= sasa <= 80.0:
        pocket_bonus = 30.0 # Es accesible y cóncavo
    elif sasa < 5.0:
        pocket_bonus = -100.0 # ENTERRADO (Como el residuo 513 antiguo) - PENALIZAR
    elif sasa > 100.0:
        pocket_bonus = -10.0 # Demasiado expuesto/plano
        
    # Fórmula Maestra
    final_score = (cons_score * 3.0) + (cent_score * 50.0) + (rigidity * 1.5) + pocket_bonus
    
    data.append({
        "Residue": res.get_resname(),
        "ID": res_id,
        "Ghost_Score": round(final_score, 3),
        "Conservation": round(cons_score, 2),
        "SASA": round(sasa, 1),
        "Rigidity": round(rigidity, 1)
    })

# Ranking
df = pd.DataFrame(data).sort_values(by="Ghost_Score", ascending=False)
df.to_csv(output_csv, index=False)
confirm_file(output_csv, "ranking de candidatos")

winner = df.iloc[0]
log_info(f"🏆 GANADOR CIENTÍFICO: {winner['Residue']}{winner['ID']}")
log_info(f"Score: {winner['Ghost_Score']} (SASA: {winner['SASA']})")
