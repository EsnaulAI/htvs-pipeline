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
import subprocess
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
    log_info,
    log_warn,
    require_file,
)

homologs_input = snakemake.input.homologs
orig_fasta = snakemake.input.original_fasta
msa_output = snakemake.output.fasta_msa
threads = snakemake.threads

require_file(homologs_input, "FASTA de homologos")
require_file(orig_fasta, "FASTA original")
ensure_parent_dir(msa_output)

log_info("📖 Leyendo homologos y preparando alineamiento...")

hits = list(SeqIO.parse(homologs_input, "fasta"))
log_info(f"Se encontraron {len(hits)} secuencias homólogas válidas.")

if len(hits) < 5:
    log_warn("⚠️ ADVERTENCIA CRÍTICA: Muy pocas secuencias para análisis evolutivo.")

log_info("🧬 Ejecutando MAFFT (Alineamiento Múltiple)...")
cmd = f"mafft --auto --quiet --thread {threads} {homologs_input} > {msa_output}"
subprocess.run(cmd, shell=True, check=True)

print("✅ Alineamiento completado.")
confirm_file(msa_output, "MSA generado")
log_info("✅ Alineamiento completado.")
