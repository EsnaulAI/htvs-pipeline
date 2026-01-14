# --- MOCK PARA DESARROLLO (Pylance no se quejará) ---
if "snakemake" not in globals():
    # Esto simula la variable snakemake para que el editor no marque error
    # y tú sepas qué estructura tiene. Solo sirve mientras editas.
    from types import SimpleNamespace
    snakemake = SimpleNamespace(
        input=SimpleNamespace(pdb="", fasta="", xml="", original_fasta="", msa=""),
        output=SimpleNamespace(
            pdb="",
            fasta="",
            xml="",
            fasta_msa="",
            pdb_conserved="",
            fasta_homologs="results/module1/unaligned_homologs.fasta",
            method_log="results/module1/homolog_method.log",
        ),
        params=SimpleNamespace(pdb_id="7KGY", chain="B", n_hits=10, e_val=0.001),
        wildcards=SimpleNamespace()
    )
# ----------------------------------------------------
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys
from pathlib import Path
import yaml
from module1_run_mmseqs import run_mmseqs_search
from logging_utils import (
    confirm_file,
    ensure_parent_dir,
    log_error,
    log_info,
    log_warn,
    require_file,
)

def load_config():
    config_path = Path(__file__).resolve().parents[2] / "config" / "config.yaml"
    with open(config_path, "r") as f:
        return yaml.safe_load(f)

def write_method_log(path, method):
    if not path:
        return
    ensure_parent_dir(path)
    with open(path, "w") as log_handle:
        log_handle.write(f"{method}\n")

def blast_xml_to_homologs(xml_input, orig_fasta, homologs_output):
    require_file(xml_input, "XML de BLAST")
    require_file(orig_fasta, "FASTA original")
    ensure_parent_dir(homologs_output)

    log_info("📖 Leyendo resultados XML y preparando homologos...")

    hits = []
    original_rec = SeqIO.read(orig_fasta, "fasta")
    hits.append(original_rec)

    try:
        with open(xml_input) as result_handle:
            blast_record = NCBIXML.read(result_handle)

            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    if hsp.align_length < len(original_rec.seq) * 0.6:
                        continue

                    rec = SeqRecord(
                        Seq(hsp.sbjct),
                        id=alignment.hit_id,
                        description=alignment.hit_def,
                    )
                    hits.append(rec)
    except Exception as e:
        log_error(f"❌ Error leyendo XML: {e}")
        sys.exit(1)

    log_info(f"Se encontraron {len(hits)} secuencias homólogas válidas.")

    SeqIO.write(hits, homologs_output, "fasta")
    confirm_file(homologs_output, "FASTA de homologos")

# Inputs de Snakemake
fasta_input = snakemake.input.fasta
xml_output = snakemake.output.xml
homologs_output = snakemake.output.fasta_homologs
method_log = snakemake.output.method_log
e_val = snakemake.params.e_val
n_hits = snakemake.params.n_hits
# Default a 'nr' si no se especifica, para máxima potencia en PyDCA
db_name = snakemake.params.get("db", "nr") 

require_file(fasta_input, "FASTA de entrada")
ensure_parent_dir(xml_output)
ensure_parent_dir(homologs_output)

log_info(f"🚀 Iniciando BLAST Remoto (NCBI) contra '{db_name}' para {fasta_input}...")
log_info(f"🎯 Objetivo: {n_hits} secuencias (E-value: {e_val})")
log_info("⏳ Esto puede tardar 10-20 minutos. No cierres la terminal...")

try:
    record = SeqIO.read(fasta_input, "fasta")

    # Usamos 'nr' (Non-Redundant) para obtener miles de secuencias y alimentar PyDCA
    result_handle = NCBIWWW.qblast(
        "blastp",
        db_name,
        record.seq,
        hitlist_size=n_hits,
        expect=e_val,
    )

    with open(xml_output, "w") as out_handle:
        out_handle.write(result_handle.read())

    result_handle.close()
    confirm_file(xml_output, "XML de resultados BLAST")
    log_info(f"✅ Resultados BLAST guardados en {xml_output}")

    blast_xml_to_homologs(xml_output, fasta_input, homologs_output)
    write_method_log(method_log, "blast")

except Exception as e:
    log_error(f"❌ Error conectando a BLAST: {e}")
    log_warn("↪️ Ejecutando MMseqs2 local como respaldo.")
    config = load_config()
    threads = getattr(snakemake, "threads", 4)
    run_mmseqs_search(fasta_input, homologs_output, config, threads, method_log_path=method_log)
    with open(xml_output, "w") as out_handle:
        out_handle.write(f"BLAST falló; se usó MMseqs2. Error: {e}\n")
    confirm_file(xml_output, "XML de resultados BLAST (fallback)")
