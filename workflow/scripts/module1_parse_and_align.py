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
        input=SimpleNamespace(pdb="", fasta="", xml="", original_fasta="", msa=""),
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
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys
import yaml

def load_config():
    config_path = Path(__file__).resolve().parents[2] / "config" / "config.yaml"
    with open(config_path, "r") as f:
        return yaml.safe_load(f)

xml_input = snakemake.input.xml
orig_fasta = snakemake.input.original_fasta
msa_output = snakemake.output.fasta_msa
threads = snakemake.threads
config = load_config()

print(f"📖 Leyendo resultados XML y preparando alineamiento...")

hits = []
# Leemos nuestra secuencia original primero
original_rec = SeqIO.read(orig_fasta, "fasta")
hits.append(original_rec)

try:
    with open(xml_input) as result_handle:
        blast_record = NCBIXML.read(result_handle)
        
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                # Filtro de Calidad:
                # El alineamiento debe cubrir al menos el 60% de nuestra proteína
                if hsp.align_length < len(original_rec.seq) * 0.6:
                    continue
                
                # Reconstruimos la secuencia
                seq = Seq(hsp.sbjct)
                rec = SeqRecord(seq, id=alignment.hit_id, description=alignment.hit_def)
                hits.append(rec)
except Exception as e:
    print(f"❌ Error leyendo XML: {e}")
    sys.exit(1)

print(f"   -> Se encontraron {len(hits)} secuencias homólogas válidas.")

if len(hits) < 5:
    print("⚠️ ADVERTENCIA CRÍTICA: Muy pocas secuencias para análisis evolutivo.")

# Guardar temporalmente
temp_fasta = config["evolution"]["unaligned_homologs"]
SeqIO.write(hits, temp_fasta, "fasta")

print(f"🧬 Ejecutando MAFFT (Alineamiento Múltiple)...")
cmd = f"mafft --auto --quiet --thread {threads} {temp_fasta} > {msa_output}"
subprocess.run(cmd, shell=True, check=True)

print("✅ Alineamiento completado.")
