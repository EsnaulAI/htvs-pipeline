import sys
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import yaml
from Bio import AlignIO
from Bio.PDB import PDBParser

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.append(str(SCRIPT_DIR))

from module1_map_entropy import map_pdb_to_msa

REPO_ROOT = SCRIPT_DIR.parents[1]
CONFIG_FILE = REPO_ROOT / "config" / "config.yaml"

config = yaml.safe_load(CONFIG_FILE.read_text())
msa_file = config["evolution"]["msa_file"]
pdb_file = config["structure"]["clean_pdb"]
chain_id = config["structure"].get("chain_id")
target_residue = config["structure"]["target_residue"]

print("🧬 Iniciando Análisis de Co-evolución (Mutual Information)...")

# 1. Cargar MSA
alignment = AlignIO.read(msa_file, "fasta")
num_seqs = len(alignment)
aln_len = alignment.get_alignment_length()

# Mapear residuo objetivo PDB -> columna MSA
parser = PDBParser(QUIET=True)
structure = parser.get_structure("Target", pdb_file)
mapping, _, chain_id_used = map_pdb_to_msa(alignment, structure, chain_id=chain_id)
if target_residue not in mapping:
    print(
        "⚠️ Advertencia: no hay correspondencia para el residuo objetivo "
        f"{target_residue} en la cadena {chain_id_used}. "
        "Es posible que esté ausente por gaps en el MSA."
    )
    sys.exit(1)

target_res_index = mapping[target_residue]
print(
    "   Residuo objetivo PDB "
    f"{target_residue} -> columna MSA {target_res_index + 1}."
)

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
target_col = msa_matrix[:, target_res_index]

# Entropía del target (para normalizar)
h_x = calc_entropy(target_col)

print("   Calculando correlaciones (esto toma unos segundos)...")
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

print("\n" + "=" * 40)
print(f"🔗 REPORTE DE CO-EVOLUCIÓN (Socios de residuo {target_residue})")
print("=" * 40)
print(f"Top 5 residuos co-evolucionando con el Target:")
for idx, score in zip(top_indices, top_scores):
    # Aquí el índice es del MSA, correspondería mapearlo al PDB en un caso ideal
    print(f"   Columna MSA {idx + 1}: NMI = {score:.4f}")

if max(top_scores) > 0.1:
    print("\n✅ CONCLUSIÓN: Existen señales de co-evolución.")
    print("   El residuo 'habla' con otros sitios a distancia.")
else:
    print("\n⚠️ CONCLUSIÓN: Baja señal de co-evolución.")
    print("   El residuo podría ser independiente o muy conservado (si no cambia, no co-evoluciona).")

# Plot rápido
plt.figure(figsize=(10, 4))
plt.plot(mi_scores)
plt.title(f"Perfil de Co-evolución para Residuo {target_residue}")
plt.xlabel("Residuo (Columna MSA)")
plt.ylabel("Información Mutua Normalizada")
plt.savefig("results/module1/coevolution_profile.png")
print("   -> Gráfica guardada: results/module1/coevolution_profile.png")
