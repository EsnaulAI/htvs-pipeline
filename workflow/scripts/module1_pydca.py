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
            msa=config["evolution"]["msa_file"],
            pdb_msa_map=config["evolution"]["pdb_msa_map"],
        ),
        output=SimpleNamespace(plot=config["analysis"]["dca_distribution"], report=config["analysis"]["dca_top_contacts"]),
        params=SimpleNamespace(
            pdb_id=config["structure"]["pdb_id"],
            chain=config["structure"]["chain_id"],
            n_hits=config["evolution"]["n_homologs"],
            e_val=config["evolution"]["e_value"],
            target_res=config["structure"]["target_residue"],
            target_res_name=config["structure"]["target_residue_name"],
        ),
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
from pathlib import Path
import json
import yaml

def load_config():
    config_path = Path(__file__).resolve().parents[2] / "config" / "config.yaml"
    with open(config_path, "r") as f:
        return yaml.safe_load(f)
from logging_utils import confirm_file, ensure_parent_dir, log_error, log_info, log_warn, require_file

msa_file = snakemake.input.msa
map_file = snakemake.input.pdb_msa_map
output_plot = snakemake.output.plot 
output_report = snakemake.output.report 

# --- CAMBIO CRÍTICO: Recibir variable desde Config ---
target_res_pdb = int(snakemake.params.target_res)
# ---------------------------------------------------
config = load_config()
target_res_name = snakemake.params.target_res_name if hasattr(snakemake.params, "target_res_name") else config["structure"]["target_residue_name"]

require_file(msa_file, "MSA de entrada")
require_file(map_file, "mapa PDB->MSA")
ensure_parent_dir(output_plot)
ensure_parent_dir(output_report)

with open(map_file, "r", encoding="utf-8") as handle:
    map_payload = json.load(handle)
mapping = {int(k): v for k, v in map_payload.get("mapping", {}).items()}
index_base = map_payload.get("index_base", 0)
chain_id_used = map_payload.get("chain_id")
if target_res_pdb not in mapping:
    log_error(
        "❌ Error: el residuo objetivo "
        f"{target_res_pdb} en la cadena {chain_id_used} no está mapeado. "
        "Revisa el mapa PDB→MSA o la presencia de gaps en el MSA."
    )
    sys.exit(1)

target_res_msa = mapping[target_res_pdb] - index_base
target_res_msa_pos = target_res_msa + 1
log_info(f"Índice base del mapa PDB→MSA: {index_base} (0-based o 1-based).")

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
            
            if abs(res1 - target_res_msa_pos) < 5 or abs(res2 - target_res_msa_pos) < 5:
                found_target = True

        if found_target:
            print(
                "\n✅ ¡BINGO! El entorno de "
                f"{target_res_name} {target_res_pdb} "
                f"(MSA {target_res_msa_pos}) aparece en los acoplamientos fuertes."
            )
        else:
            print(
                "\nℹ️ "
                f"{target_res_name} {target_res_pdb} "
                f"(MSA {target_res_msa_pos}) no está en el Top 20 global."
            )

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
    print(f"❌ Error crítico en PyDCA: {e}")
    log_error(f"❌ Error crítico en PyDCA: {e}")
    sys.exit(1)
