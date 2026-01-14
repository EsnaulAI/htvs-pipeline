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
import sys
from Bio.PDB import PDBList, PDBParser, Select, PDBIO, Polypeptide
from logging_utils import (
    confirm_file,
    ensure_parent_dir,
    log_error,
    log_info,
)

# Snakemake injects inputs/outputs/params automatically
pdb_id = snakemake.params.pdb_id
chain_target = snakemake.params.chain
out_pdb = snakemake.output.pdb
out_fasta = snakemake.output.fasta

ensure_parent_dir(out_pdb)
ensure_parent_dir(out_fasta)

log_info(f"🔬 Iniciando descarga y limpieza de {pdb_id} cadena {chain_target}...")

# 1. Descargar
pdbl = PDBList()
pdb_file = pdbl.retrieve_pdb_file(pdb_id, pdir="data/raw", file_format="pdb")

# 2. Parsear y Limpiar
class ChainSelect(Select):
    def accept_chain(self, chain):
        return chain.get_id() == chain_target
    def accept_residue(self, residue):
        # Eliminar aguas (HOH) y heterátomos que no sean parte de la proteína
        return residue.id[0] == " "

parser = PDBParser(QUIET=True)
structure = parser.get_structure(pdb_id, pdb_file)

io = PDBIO()
io.set_structure(structure)
io.save(out_pdb, select=ChainSelect())

# 3. Extraer Secuencia FASTA (FIX: Concatenar fragmentos)
ppb = Polypeptide.PPBuilder()
peptides = ppb.build_peptides(structure[0][chain_target])

if len(peptides) == 0:
    log_error("❌ Error Crítico: No se encontraron péptidos en la cadena especificada.")
    sys.exit(1)

# Aquí está el arreglo: Unimos todos los fragmentos en una sola string
full_seq = "".join([str(pp.get_sequence()) for pp in peptides])

log_info(f"ℹ️  Nota: La estructura tiene {len(peptides)} fragmentos. Se han unido para el análisis BLAST.")

with open(out_fasta, "w") as f:
    f.write(f">{pdb_id}_{chain_target}\n{full_seq}\n")

confirm_file(out_pdb, "PDB limpio")
confirm_file(out_fasta, "FASTA de salida")
log_info("✅ Estructura limpia y secuencia unificada extraída.")
