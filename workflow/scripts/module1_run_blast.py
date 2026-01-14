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
from Bio.Blast import NCBIWWW
from Bio import SeqIO
import sys
from logging_utils import (
    confirm_file,
    ensure_parent_dir,
    log_error,
    log_info,
    require_file,
)

# Inputs de Snakemake
fasta_input = snakemake.input.fasta
xml_output = snakemake.output.xml
e_val = snakemake.params.e_val
n_hits = snakemake.params.n_hits
# Default a 'nr' si no se especifica, para máxima potencia en PyDCA
db_name = snakemake.params.get("db", "nr") 

require_file(fasta_input, "FASTA de entrada")
ensure_parent_dir(xml_output)

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
        expect=e_val
    )

    with open(xml_output, "w") as out_handle:
        out_handle.write(result_handle.read())

    result_handle.close()
    confirm_file(xml_output, "XML de resultados BLAST")
    log_info(f"✅ Resultados BLAST guardados en {xml_output}")

except Exception as e:
    log_error(f"❌ Error conectando a BLAST: {e}")
    sys.exit(1)
