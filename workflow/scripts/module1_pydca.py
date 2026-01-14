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

from pydca.meanfield_dca import meanfield_dca
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from logging_utils import (
    confirm_file,
    ensure_parent_dir,
    log_error,
    log_info,
    log_warn,
    require_file,
)

msa_file = snakemake.input.msa
output_plot = snakemake.output.plot 
output_report = snakemake.output.report 

# --- CAMBIO CRÍTICO: Recibir variable desde Config ---
target_res_pdb = int(snakemake.params.target_res)
# ---------------------------------------------------

require_file(msa_file, "MSA de entrada")
ensure_parent_dir(output_plot)
ensure_parent_dir(output_report)

log_info("🧠 Iniciando PyDCA (Direct Coupling Analysis) - Campo Medio...")

# 0. VERIFICACIÓN DE SEGURIDAD (CRÍTICO)
records = list(SeqIO.parse(msa_file, "fasta"))
num_seqs = len(records)
log_info(f"📊 Secuencias detectadas para análisis: {num_seqs}")

if num_seqs < 5:
    log_warn("⚠️  ADVERTENCIA: Insuficientes secuencias para análisis evolutivo (DCA).")
    log_warn("Se requieren homólogos variados para detectar co-evolución.")
    log_warn("-> Generando archivos vacíos para completar el pipeline sin errores.")
    
    with open(output_report, "w") as f:
        f.write("Error,Residue1,Residue2,Score\n")
        f.write(f"SKIPPED,0,0,0.0\n")
        f.write(f"# Reason: Only {num_seqs} sequences found. Need diversity for DCA.\n")
    
    plt.figure()
    plt.text(0.1, 0.5, f"DCA Omitido: Solo {num_seqs} secuencia(s).\nSe requiere mas diversidad.", fontsize=12)
    plt.savefig(output_plot)
    
    confirm_file(output_report, "reporte DCA")
    confirm_file(output_plot, "gráfica DCA")
    sys.exit(0)

# 1. Instanciar DCA
try:
    mfdca = meanfield_dca.MeanFieldDCA(
        msa_file,
        biomolecule="protein", 
        pseudocount=0.5,
        seqid=0.8
    )

    # 2. Calcular Matriz de Acoplamiento
    log_info("Calculando matriz inversa (esto consume CPU)...")
    fnapc = mfdca.compute_sorted_FN_APC() 

    # 3. Analizar Resultados
    log_info("Analizando top interacciones fuertes...")
    top_interactions = fnapc[:20] 

    with open(output_report, "w") as f:
        f.write("Residue1,Residue2,Score_DCA\n")
        found_target = False
        
        log_info("🏆 Top Conexiones Evolutivas Reales (DCA):")
        for pair in top_interactions:
            # pair es ((i, chain), (j, chain), score)
            res1, res2, score = pair[0][0], pair[1][0], pair[2]
            f.write(f"{res1},{res2},{score}\n")
            log_info(f"Res {res1} <--> Res {res2} : Score {score:.4f}")
            
            if abs(res1 - target_res_pdb) < 5 or abs(res2 - target_res_pdb) < 5:
                found_target = True

        if found_target:
            log_info(f"✅ ¡BINGO! El entorno de GLU {target_res_pdb} aparece en los acoplamientos fuertes.")
        else:
            log_info(f"ℹ️ GLU {target_res_pdb} no está en el Top 20 global.")

    # 4. Plot de Contactos
    data = np.array([x[2] for x in fnapc])
    plt.figure()
    plt.plot(data)
    plt.title("Distribución de Scores DCA (Ley de Potencia)")
    plt.xlabel("Ranking de Pares")
    plt.ylabel("Fuerza de Acoplamiento")
    plt.savefig(output_plot)
    confirm_file(output_report, "reporte DCA")
    confirm_file(output_plot, "gráfica DCA")
    log_info(f"-> Gráfica guardada: {output_plot}")

except Exception as e:
    log_error(f"❌ Error crítico en PyDCA: {e}")
    sys.exit(1)
