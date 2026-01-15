# --- MOCK PARA DESARROLLO (Pylance no se quejará) ---
if "snakemake" not in globals():
    # Esto simula la variable snakemake para que el editor no marque error
    # y tú sepas qué estructura tiene. Solo sirve mientras editas.
    from pathlib import Path
    from types import SimpleNamespace
    import yaml

    def load_config():
        config_path = Path(__file__).resolve().parents[2] / "config" / "config.yaml"
        with open(config_path, "r") as f:
            return yaml.safe_load(f)

    config = load_config()
    snakemake = SimpleNamespace(
        input=SimpleNamespace(pdb="", fasta="", xml="", original_fasta="", msa="", homologs=config["evolution"]["unaligned_homologs"]),
        output=SimpleNamespace(pdb="", fasta="", xml="", fasta_msa=config["evolution"]["msa_file"], pdb_conserved=""),
        params=SimpleNamespace(
            pdb_id=config["structure"]["pdb_id"],
            chain=config["structure"]["chain_id"],
            n_hits=config["evolution"]["n_homologs"],
            e_val=config["evolution"]["e_value"],
        ),
        wildcards=SimpleNamespace()
    )
# ----------------------------------------------------
import os
import subprocess
import sys
import tempfile
from pathlib import Path
from Bio import SeqIO
import yaml

def load_config():
    config_path = Path(__file__).resolve().parents[2] / "config" / "config.yaml"
    with open(config_path, "r") as f:
        return yaml.safe_load(f)
from logging_utils import (
    confirm_file,
    ensure_parent_dir,
    log_error,
    log_info,
    log_warn,
    require_file,
)

homologs_input = snakemake.input.homologs
orig_fasta = snakemake.input.original_fasta
msa_output = snakemake.output.fasta_msa
threads = snakemake.threads
config = load_config()
min_msa_sequences = config["evolution"].get("min_msa_sequences", 5)
min_msa_length = config["evolution"].get("min_msa_length", 50)

require_file(homologs_input, "FASTA de homologos")
require_file(orig_fasta, "FASTA original")
ensure_parent_dir(msa_output)

log_info("📖 Leyendo homologos y preparando alineamiento...")

hits = list(SeqIO.parse(homologs_input, "fasta"))
log_info(f"Se encontraron {len(hits)} secuencias homólogas válidas.")

if len(hits) < 5:
    log_warn("⚠️ ADVERTENCIA CRÍTICA: Muy pocas secuencias para análisis evolutivo.")

log_info("🧬 Ejecutando MAFFT (Alineamiento Múltiple)...")
temp_fasta = homologs_input
temp_log_path = None
try:
    with tempfile.NamedTemporaryFile(mode="w", suffix=".log", delete=False) as log_handle:
        temp_log_path = log_handle.name
        cmd = ["mafft", "--auto", "--quiet", "--thread", str(threads), temp_fasta]
        with open(msa_output, "w") as msa_handle:
            subprocess.run(cmd, check=True, stdout=msa_handle, stderr=log_handle)
except subprocess.CalledProcessError as exc:
    log_error(f"❌ Error ejecutando MAFFT: {exc}")
    sys.exit(1)
finally:
    if temp_log_path:
        os.unlink(temp_log_path)

print("✅ Alineamiento completado.")
confirm_file(msa_output, "MSA generado")
log_info("✅ Alineamiento completado.")

alignment_records = list(SeqIO.parse(msa_output, "fasta"))
if len(alignment_records) < min_msa_sequences:
    log_error(
        "❌ MSA insuficiente: "
        f"{len(alignment_records)} secuencias (mínimo {min_msa_sequences})."
    )
    sys.exit(1)

alignment_length = len(alignment_records[0].seq) if alignment_records else 0
if alignment_length <= min_msa_length:
    log_error(
        "❌ MSA insuficiente: "
        f"longitud {alignment_length} (debe ser > {min_msa_length})."
    )
    sys.exit(1)
