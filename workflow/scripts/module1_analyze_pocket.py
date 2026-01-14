# --- MOCK PARA DESARROLLO (Pylance no se quejará) ---
if "snakemake" not in globals():
    from types import SimpleNamespace
    snakemake = SimpleNamespace(
        input=SimpleNamespace(msa="results/module1/alignment.fasta"),
        output=SimpleNamespace(plot="results/module1/dca_distribution.png", report="results/module1/dca_top_contacts.csv"),
        params=SimpleNamespace(pdb_id="7KGY", chain="B", n_hits=10, e_val=0.001, target_res=513),
        wildcards=SimpleNamespace(),
        threads=4
    )
# ----------------------------------------------------
# workflow/scripts/module1_analyze_pocket.py

from Bio.PDB import PDBParser, NeighborSearch, Selection
import numpy as np
import sys
from logging_utils import (
    confirm_file,
    ensure_parent_dir,
    log_error,
    log_info,
    require_file,
)

# Inputs
pdb_file = snakemake.input.pdb
chain_id = snakemake.params.chain
target_res_num = int(snakemake.params.target_res)
output_report = snakemake.output.report

require_file(pdb_file, "PDB conservado")
ensure_parent_dir(output_report)

log_info(f"🏘️ Analizando el vecindario evolutivo alrededor de {chain_id}:{target_res_num}...")

parser = PDBParser(QUIET=True)
structure = parser.get_structure("Target", pdb_file)

# 1. Obtener átomos de la cadena correcta
atom_list = [atom for atom in structure[0][chain_id].get_atoms()]
ns = NeighborSearch(atom_list)

# 2. Encontrar el residuo central
target_residue = None
for res in structure[0][chain_id]:
    if res.id[1] == target_res_num:
        target_residue = res
        break

if not target_residue:
    log_error("❌ Error: Target no encontrado.")
    sys.exit(1)

# Usamos el Carbono Alpha o un átomo central para buscar vecinos
center_atom = target_residue['CA'] if 'CA' in target_residue else target_residue.child_list[0]

# 3. Buscar vecinos a 6 Angstroms (Radio de contacto de un fármaco)
neighbors = ns.search(center_atom.get_coord(), 6.0, level='R')

# 4. Analizar Conservación (Extraída del B-factor)
pocket_scores = []
log_info("📋 RESIDUOS DEL BOLSILLO:")
log_info(f"{'Residuo':<10} {'Conservación (0-4.3)':<20} {'Estado'}")
log_info("-" * 45)

with open(output_report, "w") as f:
    f.write("Residue,Conservation_Score,Status\n")
    
    for res in neighbors:
        # El score se guardó en el b_factor en el paso map_conservation
        # Tomamos el promedio de los átomos del residuo (debería ser igual para todos)
        score = res['CA'].get_bfactor() if 'CA' in res else list(res.get_atoms())[0].get_bfactor()
        pocket_scores.append(score)
        
        status = "Variable ⚠️"
        if score > 2.0: status = "Estable"
        if score > 3.5: status = "Crítico 🔥"
        
        line = f"{res.get_resname()}{res.id[1]:<5} {score:<20.2f} {status}"
        log_info(line)
        f.write(f"{res.get_resname()}{res.id[1]},{score:.2f},{status}\n")

# 5. Veredicto Final
avg_pocket = np.mean(pocket_scores)
log_info("-" * 45)
log_info(f"📊 Puntuación Global del Bolsillo: {avg_pocket:.2f} / 4.32")

conclusion = ""
if avg_pocket > 3.0:
    conclusion = "✅ BOLSILLO DE ALTA CONFIANZA: Ideal para fármacos de amplio espectro."
elif avg_pocket > 1.5:
    conclusion = "⚠️ BOLSILLO HÍBRIDO: El fármaco debe anclarse a los residuos conservados."
else:
    conclusion = "⛔ BOLSILLO HIPER-VARIABLE: Alto riesgo de resistencia."

log_info(f"-> {conclusion}")
with open(output_report, "a") as f:
    f.write(f"\n# GLOBAL SCORE: {avg_pocket:.2f}\n")
    f.write(f"# CONCLUSION: {conclusion}\n")

confirm_file(output_report, "reporte del bolsillo")
