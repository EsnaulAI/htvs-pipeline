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
# workflow/scripts/module1_check_druggability.py

from Bio.PDB import PDBParser, SASA
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

require_file(pdb_file, "PDB de entrada")
ensure_parent_dir(output_report)

log_info(f"💧 Calculando Farmacabilidad (SASA) para {chain_id}:{target_res_num}...")

try:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("Target", pdb_file)
    target_chain = structure[0][chain_id]

    # 1. Calcular SASA (Solvent Accessible Surface Area)
    # Requiere que tengas 'freesasa' o usar el modulo interno de Biopython
    sr = SASA.ShrakeRupley()
    sr.compute(structure, level="R") # Nivel Residuo

    target_sasa = 0
    all_sasa = []

    for residue in target_chain:
        # Biopython agrega el atributo .sasa al residuo
        if hasattr(residue, 'sasa'):
            sasa = residue.sasa
            all_sasa.append(sasa)
            if residue.id[1] == target_res_num:
                target_sasa = sasa

    mean_sasa = np.mean(all_sasa) if all_sasa else 0

    # 2. Diagnóstico
    log_info("📊 DIAGNÓSTICO DE ACCESIBILIDAD:")
    log_info(f"SASA del Target: {target_sasa:.2f} Å²")
    log_info(f"Promedio Proteína: {mean_sasa:.2f} Å²")

    conclusion = ""
    if target_sasa < 5:
        conclusion = "⛔ BURIED (Enterrado): Sitio profundo. Requiere dinámica para abrirse."
    elif 5 <= target_sasa <= 30:
        conclusion = "✅ POCKET (Bolsillo): Accesibilidad ideal para sitio críptico (Cryptic Pocket)."
    else:
        conclusion = "⚠️ EXPOSED (Expuesto): Superficie muy abierta. Difícil selectividad."

    log_info(f"-> {conclusion}")

    with open(output_report, "w") as f:
        f.write(f"Target Residue: {target_res_num}\n")
        f.write(f"SASA_Score: {target_sasa:.2f}\n")
        f.write(f"Global_Mean_SASA: {mean_sasa:.2f}\n")
        f.write(f"Conclusion: {conclusion}\n")
        f.write("Ghost_Score_Pocket_Weights:\n")
        f.write("  Pocket_Member_Weight: 20.0\n")
        f.write("  SASA_Weight: 10.0\n")
        f.write(
            "Justification: Se priorizan residuos confirmados por fpocket como parte de un bolsillo "
            "real, mientras que la SASA normalizada mantiene el criterio de accesibilidad.\n"
        )

    confirm_file(output_report, "reporte de farmacabilidad")
except Exception as e:
    log_error(f"❌ Error calculando SASA: {e}")
    # Generar reporte de error para no romper el pipeline
    with open(output_report, "w") as f:
        f.write("Error calculation SASA\n")
    confirm_file(output_report, "reporte de farmacabilidad")
