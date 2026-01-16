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
        input=SimpleNamespace(msa=config["evolution"]["msa_file"]),
        output=SimpleNamespace(report=config["evolution"]["msa_qc_report"]),
        params=SimpleNamespace(),
        wildcards=SimpleNamespace(),
    )
# ----------------------------------------------------

from itertools import combinations
from Bio import AlignIO

from logging_utils import (
    confirm_file,
    ensure_parent_dir,
    log_info,
    log_error,
    log_warn,
    require_file,
)


def calculate_average_coverage(alignment):
    aln_len = alignment.get_alignment_length()
    if aln_len == 0:
        return 0.0
    coverages = []
    for record in alignment:
        seq = str(record.seq)
        non_gaps = sum(1 for aa in seq if aa != "-")
        coverages.append(non_gaps / aln_len)
    return sum(coverages) / len(coverages) if coverages else 0.0


def calculate_gap_percentages(alignment):
    num_seqs = len(alignment)
    aln_len = alignment.get_alignment_length()
    gap_percentages = []
    for col_idx in range(aln_len):
        column = alignment[:, col_idx]
        gaps = sum(1 for aa in column if aa == "-")
        gap_percent = (gaps / num_seqs) * 100 if num_seqs else 0.0
        gap_percentages.append(gap_percent)
    return gap_percentages


def calculate_mean_pairwise_identity(alignment):
    sequences = [str(record.seq) for record in alignment]
    if len(sequences) < 2:
        return 0.0

    identities = []
    for seq_a, seq_b in combinations(sequences, 2):
        matches = 0
        comparable = 0
        for aa, bb in zip(seq_a, seq_b):
            if aa == "-" or bb == "-":
                continue
            comparable += 1
            if aa == bb:
                matches += 1
        if comparable:
            identities.append(matches / comparable)

    return sum(identities) / len(identities) if identities else 0.0


msa_file = snakemake.input.msa
output_report = snakemake.output.report
min_sequences = 50
min_coverage = 0.6

require_file(msa_file, "MSA de entrada")
ensure_parent_dir(output_report)

log_info("🔎 Ejecutando QC del MSA...")

alignment = AlignIO.read(msa_file, "fasta")
num_seqs = len(alignment)
aln_len = alignment.get_alignment_length()

if num_seqs == 0 or aln_len == 0:
    log_warn("⚠️ MSA vacío o sin columnas. Generando reporte mínimo.")
    with open(output_report, "w") as handle:
        handle.write("metric,value\n")
        handle.write(f"num_seqs,{num_seqs}\n")
        handle.write(f"alignment_length,{aln_len}\n")
    confirm_file(output_report, "reporte QC MSA")
    log_error(
        f"El MSA no cumple con los mínimos requeridos (secuencias >= {min_sequences}, "
        f"cobertura >= {min_coverage:.2f}). Se detectó num_seqs={num_seqs}, "
        f"alignment_length={aln_len}."
    )
    raise SystemExit(1)
else:
    avg_coverage = calculate_average_coverage(alignment)
    gap_percentages = calculate_gap_percentages(alignment)
    mean_identity = calculate_mean_pairwise_identity(alignment)

    if num_seqs < min_sequences or avg_coverage < min_coverage:
        log_error(
            f"El MSA no cumple con los mínimos requeridos (secuencias >= {min_sequences}, "
            f"cobertura >= {min_coverage:.2f}). Se detectó num_seqs={num_seqs}, "
            f"average_coverage={avg_coverage:.4f}."
        )
        raise SystemExit(1)

    with open(output_report, "w") as handle:
        handle.write("metric,value\n")
        handle.write(f"num_seqs,{num_seqs}\n")
        handle.write(f"alignment_length,{aln_len}\n")
        handle.write(f"average_coverage,{avg_coverage:.4f}\n")
        handle.write(f"mean_identity,{mean_identity:.4f}\n")
        handle.write("\ncolumn,gap_percent\n")
        for idx, gap_percent in enumerate(gap_percentages, start=1):
            handle.write(f"{idx},{gap_percent:.2f}\n")

    confirm_file(output_report, "reporte QC MSA")
    log_info(f"✅ Reporte QC generado: {output_report}")
