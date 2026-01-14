# --- MOCK PARA DESARROLLO ---
if "snakemake" not in globals():
    from types import SimpleNamespace
    snakemake = SimpleNamespace(
        input=SimpleNamespace(pdb="results/module1/target_conserved.pdb"),
        output=SimpleNamespace(report="results/module1/network_report.txt"),
        params=SimpleNamespace(chain="B", target_res=513)
    )
# ----------------------------

import networkx as nx
import numpy as np
from Bio.PDB import PDBParser
import sys
from logging_utils import (
    confirm_file,
    ensure_parent_dir,
    log_info,
    require_file,
)

# Inputs
pdb_file = snakemake.input.pdb
chain_id = snakemake.params.chain

# --- CAMBIO CRÍTICO: Recibir variable desde Config ---
target_residue = int(snakemake.params.target_res)
# ---------------------------------------------------

output_report = snakemake.output.report

require_file(pdb_file, "PDB de entrada")
ensure_parent_dir(output_report)

log_info(f"🕸️ Iniciando Análisis de Redes (Graph Theory) en cadena {chain_id}...")

parser = PDBParser(QUIET=True)
structure = parser.get_structure("Target", pdb_file)
chain = structure[0][chain_id]

# 1. Construir el Grafo (Residue Interaction Network)
G = nx.Graph()
residues = list(chain.get_residues())
residues = [r for r in residues if r.id[0] == " "] # Solo aminoácidos estándar

# Añadir nodos
for r in residues:
    G.add_node(r.id[1], resname=r.get_resname())

# Añadir aristas (Contactos < 8 Angstroms entre Carbonos Alpha)
log_info("Calculando contactos físicos...")
for i, r1 in enumerate(residues):
    for j, r2 in enumerate(residues):
        if i >= j: continue 
        
        try:
            d = r1['CA'] - r2['CA']
            if d < 8.0:
                G.add_edge(r1.id[1], r2.id[1])
        except KeyError:
            continue

log_info(f"Grafo construido: {G.number_of_nodes()} nodos, {G.number_of_edges()} conexiones.")

# 2. Calcular Centralidad
log_info("🧮 Calculando 'Betweenness Centrality'...")
centrality = nx.betweenness_centrality(G)

# 3. Analizar tu Blanco
my_score = centrality.get(target_residue, 0)
max_score = max(centrality.values()) if centrality else 0
rank = sorted(centrality.values(), reverse=True).index(my_score) + 1 if centrality else 0

# --- Generar Reporte ---
with open(output_report, "w") as f:
    f.write(f"REPORTE DE RED PARA RESIDUO {target_residue}\n")
    f.write("="*40 + "\n")
    f.write(f"Centralidad: {my_score:.4f} (Max: {max_score:.4f})\n")
    f.write(f"Ranking: #{rank} de {len(residues)}\n")
    
    log_info("=" * 40)
    log_info(f"📊 REPORTE DE RED PARA GLU {target_residue}")
    log_info(f"Centralidad: {my_score:.4f} (Max en proteína: {max_score:.4f})")
    log_info(f"Ranking: #{rank} de {len(residues)} residuos")

    if rank < len(residues) * 0.1:
        msg = "✅ CONCLUSIÓN: Es un HUB de comunicación (Top 10%). Confirmado ALOSTÉRICO."
    else:
        msg = "⚠️ CONCLUSIÓN: No es un hub central estructural."
    
    f.write(msg + "\n")
    log_info(msg)

confirm_file(output_report, "reporte de red")
