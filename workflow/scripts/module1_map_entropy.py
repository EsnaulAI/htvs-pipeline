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
import math
from Bio import AlignIO
from Bio.PDB import PDBParser, PDBIO

pdb_input = snakemake.input.pdb
msa_input = snakemake.input.msa
pdb_output = snakemake.output.pdb_conserved

def calculate_shannon_entropy(alignment):
    """Calcula la entropía de Shannon por columna del alineamiento.
    Baja entropía = Alta conservación."""
    n_seqs = len(alignment)
    aln_len = alignment.get_alignment_length()
    entropy_scores = []

    for col_idx in range(aln_len):
        column = alignment[:, col_idx]
        # Contar frecuencias de aminoácidos (ignorando gaps '-')
        counts = {}
        non_gap_count = 0
        for residue in column:
            if residue == '-' or residue == 'X':
                continue
            counts[residue] = counts.get(residue, 0) + 1
            non_gap_count += 1
        
        if non_gap_count == 0:
            entropy_scores.append(0) # Columna vacía
            continue

        entropy = 0
        for residue in counts:
            p = counts[residue] / non_gap_count
            entropy -= p * math.log2(p)
        
        # Invertir puntuación: Queremos que Alto Valor = Alta Conservación
        # La entropía máxima para 20 AA es log2(20) ≈ 4.32
        conservation = 4.32 - entropy
        if conservation < 0: conservation = 0
        entropy_scores.append(conservation)
        
    return entropy_scores

print("🧮 Calculando conservación evolutiva...")
alignment = AlignIO.read(msa_input, "fasta")
scores = calculate_shannon_entropy(alignment)

# Mapear al PDB
parser = PDBParser(QUIET=True)
structure = parser.get_structure("Target", pdb_input)
# Asumimos que la primera secuencia del MSA es nuestra PDB (porque la pusimos primera)
# Necesitamos mapear índice MSA -> Residuo PDB (saltando gaps en la seq 1)

target_seq_in_msa = alignment[0].seq
pdb_residues = list(structure.get_residues())

msa_index = 0
pdb_index = 0

print("🎨 Inyectando puntuaciones en el factor B del PDB...")

for msa_char in target_seq_in_msa:
    score = scores[msa_index]
    
    if msa_char != '-':
        # Si no es un gap en nuestra secuencia, corresponde a un residuo del PDB
        if pdb_index < len(pdb_residues):
            residue = pdb_residues[pdb_index]
            # Asignar score al B-factor de cada átomo del residuo
            for atom in residue:
                atom.set_bfactor(score)
            pdb_index += 1
            
    msa_index += 1

# Guardar PDB con datos inyectados
io = PDBIO()
io.set_structure(structure)
io.save(pdb_output)

print(f"✅ ¡Éxito! Archivo generado: {pdb_output}")
print("   -> Abre este archivo en PyMOL y colorea por B-factor para ver la conservación.")