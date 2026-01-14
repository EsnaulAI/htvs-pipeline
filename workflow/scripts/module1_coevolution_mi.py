from pathlib import Path
import sys
from collections import Counter

import matplotlib.pyplot as plt
import numpy as np
import yaml
from Bio import AlignIO

from logging_utils import confirm_file, ensure_parent_dir, log_info, log_warn, require_file


def load_config():
    config_path = Path(__file__).resolve().parents[2] / "config" / "config.yaml"
    with open(config_path, "r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def load_pdb_msa_map(map_path):
    with open(map_path, "r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle)
    mapping = {int(k): v for k, v in data.get("mapping", {}).items()}
    return mapping, data.get("chain_id"), data.get("index_base", 0)


# Inputs
msa_file = config["evolution"]["msa_file"]
target_res_name = config["structure"]["target_residue_name"]
output_plot = config["analysis"]["coevolution_profile"]
msa_file = "results/module1/alignment.fasta"
try:
    target_residue = int(snakemake.params.target_res)
except NameError:
    target_residue = config["structure"]["target_residue"]
target_res_index = target_residue - 1  # Ajuste de índice (PDB -> Array si empieza en 1, pero depende del mapeo. Usaremos aproximado)
# NOTA: En un pipeline real estricto, debemos alinear índice PDB <-> índice MSA.
# Aquí asumiremos que el MSA mantiene la numeración aproximada tras recortar gaps principales.
output_plot = "results/module1/coevolution_profile.png"
def calc_entropy(col_data):
    counts = np.bincount(col_data, minlength=21)
    probs = counts[counts > 0] / len(col_data)
    return -np.sum(probs * np.log2(probs))


def calc_mi(col_i, col_j):
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


if "snakemake" not in globals():
    from types import SimpleNamespace

    config = load_config()
    snakemake = SimpleNamespace(
        input=SimpleNamespace(
            msa=config["evolution"]["msa_file"],
            pdb_msa_map=config["evolution"]["pdb_msa_map"],
        ),
        output=SimpleNamespace(plot=config["analysis"]["coevolution_profile"]),
        params=SimpleNamespace(target_res=config["structure"]["target_residue"]),
    )
else:
    config = load_config()

msa_file = snakemake.input.msa
output_plot = snakemake.output.plot
map_file = snakemake.input.pdb_msa_map
target_residue = int(snakemake.params.target_res)
if "snakemake" in globals():
    msa_file = snakemake.input.msa
    output_plot = snakemake.output.plot
    map_file = snakemake.input.pdb_msa_map
    target_residue = int(snakemake.params.target_res)
else:
    config = load_config()
    msa_file = config["evolution"]["msa_file"]
    output_plot = config["analysis"]["coevolution_profile"]
    map_file = config["evolution"]["pdb_msa_map"]
    target_residue = config["structure"]["target_residue"]

require_file(msa_file, "MSA de entrada")
require_file(map_file, "mapa PDB->MSA")
ensure_parent_dir(output_plot)

log_info("🧬 Iniciando Análisis de Co-evolución (Mutual Information)...")

alignment = AlignIO.read(msa_file, "fasta")
num_seqs = len(alignment)
aln_len = alignment.get_alignment_length()

mapping, chain_id_used, index_base = load_pdb_msa_map(map_file)
if target_residue not in mapping:
    log_warn(
        "⚠️ Advertencia: no hay correspondencia para el residuo objetivo "
        f"{target_residue} en la cadena {chain_id_used}. "
        "Es posible que esté ausente por gaps en el MSA."
    )
    sys.exit(1)

target_res_index = mapping[target_residue] - index_base
log_info(
    "Residuo objetivo PDB "
    f"{target_residue} -> columna MSA {target_res_index + 1}."
)

aa_map = {aa: i for i, aa in enumerate("ACDEFGHIKLMNPQRSTVWY-")}
msa_matrix = np.zeros((num_seqs, aln_len), dtype=int)

for i, record in enumerate(alignment):
    seq_encoded = [aa_map.get(aa, 20) for aa in record.seq]
    msa_matrix[i, :] = seq_encoded

mi_scores = []
target_col = msa_matrix[:, target_res_index]
h_x = calc_entropy(target_col)

log_info("Calculando correlaciones (esto toma unos segundos)...")
for j in range(aln_len):
    if j == target_res_index:
        mi_scores.append(0)
        continue

    other_col = msa_matrix[:, j]
    mi = calc_mi(target_col, other_col)

    h_y = calc_entropy(other_col)
    if h_x + h_y == 0:
        nmi = 0
    else:
        nmi = (2 * mi) / (h_x + h_y)

    mi_scores.append(nmi)

top_indices = np.argsort(mi_scores)[-5:][::-1]
top_scores = [mi_scores[i] for i in top_indices]

log_info("=" * 40)
log_info(f"🔗 REPORTE DE CO-EVOLUCIÓN (Residuo {target_residue})")
log_info("=" * 40)
log_info("Top 5 residuos co-evolucionando con el Target:")
for idx, score in zip(top_indices, top_scores):
    log_info(f"Columna MSA {idx + 1}: NMI = {score:.4f}")

if max(top_scores) > 0.1:
    log_info("✅ CONCLUSIÓN: Existen señales de co-evolución.")
else:
    log_warn("⚠️ CONCLUSIÓN: Baja señal de co-evolución.")

plt.figure(figsize=(10, 4))
plt.plot(mi_scores)
plt.title(f"Perfil de Co-evolución para Residuo {target_residue}")
plt.xlabel("Residuo (Columna MSA)")
plt.ylabel("Información Mutua Normalizada")
plt.savefig(output_plot)
confirm_file(output_plot, "gráfica co-evolución")
log_info(f"-> Gráfica guardada: {output_plot}")
