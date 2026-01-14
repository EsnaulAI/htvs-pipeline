from pathlib import Path
from collections import Counter
from Bio import AlignIO
import matplotlib.pyplot as plt
import numpy as np
import sys
import yaml

def load_config():
    config_path = Path(__file__).resolve().parents[2] / "config" / "config.yaml"
    with open(config_path, "r") as f:
        return yaml.safe_load(f)

config = load_config()
from logging_utils import (
    confirm_file,
    ensure_parent_dir,
    log_info,
    log_warn,
    require_file,
)

# Inputs
msa_file = config["evolution"]["msa_file"]
target_res_index = config["structure"]["target_residue"] - 1
target_res_name = config["structure"]["target_residue_name"]
output_plot = config["analysis"]["coevolution_profile"]
# NOTA: En un pipeline real estricto, debemos alinear índice PDB <-> índice MSA. 
# Aquí asumiremos que el MSA mantiene la numeración aproximada tras recortar gaps principales.
output_plot = "results/module1/coevolution_profile.png"

require_file(msa_file, "MSA de entrada")
ensure_parent_dir(output_plot)

log_info(f"🧬 Iniciando Análisis de Co-evolución (Mutual Information) para columna ~{target_res_index}...")

# 1. Cargar MSA
alignment = AlignIO.read(msa_file, "fasta")
num_seqs = len(alignment)
aln_len = alignment.get_alignment_length()

# Convertir a matriz numpy para velocidad
# Codificamos AA como enteros
aa_map = {aa: i for i, aa in enumerate("ACDEFGHIKLMNPQRSTVWY-")}
msa_matrix = np.zeros((num_seqs, aln_len), dtype=int)

for i, record in enumerate(alignment):
    # Solo tomamos residuos estándar, otros a gap
    seq_encoded = [aa_map.get(aa, 20) for aa in record.seq]
    msa_matrix[i, :] = seq_encoded

def calc_entropy(col_data):
    counts = np.bincount(col_data, minlength=21)
    probs = counts[counts > 0] / len(col_data)
    return -np.sum(probs * np.log2(probs))

def calc_mi(col_i, col_j):
    # Probabilidad conjunta
    # Un truco rápido para pares: col_i * 100 + col_j (hashing simple)
    # Pero para velocidad en python puro, usaremos zip
    pairs = list(zip(col_i, col_j))
    counts = Counter(pairs)
    total = len(pairs)
    mi = 0.0
    for pair, count in counts.items():
        p_xy = count / total
        p_x = np.sum(col_i == pair[0]) / total
        p_y = np.sum(col_j == pair[1]) / total
        mi += p_xy * np.log2(p_xy / (p_x * p_y))
    return mi

# 2. Calcular MI para el Target contra todos
mi_scores = []
target_col = msa_matrix[:, target_res_index] # Asumimos alineación directa por ahora

# Entropía del target (para normalizar)
h_x = calc_entropy(target_col)

log_info("Calculando correlaciones (esto toma unos segundos)...")
for j in range(aln_len):
    if j == target_res_index:
        mi_scores.append(0)
        continue
        
    other_col = msa_matrix[:, j]
    mi = calc_mi(target_col, other_col)
    
    # Normalizar (NMI)
    h_y = calc_entropy(other_col)
    if h_x + h_y == 0:
        nmi = 0
    else:
        nmi = (2 * mi) / (h_x + h_y) # Media armónica simple
        
    mi_scores.append(nmi)

# 3. Análisis de Resultados
top_indices = np.argsort(mi_scores)[-5:][::-1] # Top 5
top_scores = [mi_scores[i] for i in top_indices]

print("\n" + "="*40)
print(f"🔗 REPORTE DE CO-EVOLUCIÓN (Socios de {target_res_name} {target_res_index + 1})")
print("="*40)
print(f"Top 5 residuos co-evolucionando con el Target:")
log_info("=" * 40)
log_info(f"🔗 REPORTE DE CO-EVOLUCIÓN (Socios de GLU {target_res_index + 1})")
log_info("=" * 40)
log_info("Top 5 residuos co-evolucionando con el Target:")
for idx, score in zip(top_indices, top_scores):
    # Aquí el índice es del MSA, correspondería mapearlo al PDB en un caso ideal
    log_info(f"Columna MSA {idx + 1}: NMI = {score:.4f}")

if max(top_scores) > 0.1:
    log_info("✅ CONCLUSIÓN: Existen señales de co-evolución.")
    log_info("El residuo 'habla' con otros sitios a distancia.")
else:
    log_warn("⚠️ CONCLUSIÓN: Baja señal de co-evolución.")
    log_warn("El residuo podría ser independiente o muy conservado (si no cambia, no co-evoluciona).")

# Plot rápido
plt.figure(figsize=(10, 4))
plt.plot(mi_scores)
plt.title(f"Perfil de Co-evolución para Residuo {target_res_index + 1}")
plt.xlabel("Residuo (Columna MSA)")
plt.ylabel("Información Mutua Normalizada")
plt.savefig(output_plot)
print(f"   -> Gráfica guardada: {output_plot}")
confirm_file(output_plot, "gráfica co-evolución")
log_info(f"-> Gráfica guardada: {output_plot}")
