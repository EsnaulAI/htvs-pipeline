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
from Bio import AlignIO, pairwise2
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.Polypeptide import is_aa, three_to_one

def extract_chain_sequence(structure, chain_id=None):
    model = next(structure.get_models())
    if chain_id:
        if chain_id not in model:
            raise ValueError(f"Cadena '{chain_id}' no encontrada en el PDB.")
        chain = model[chain_id]
    else:
        chain = next(model.get_chains())
    residues = [res for res in chain.get_residues() if is_aa(res, standard=True)]
    sequence = "".join(three_to_one(res.get_resname()) for res in residues)
    return residues, sequence, chain.id


def map_pdb_to_msa(alignment, structure, chain_id=None, gap_char="-"):
    """Mapea residuos del PDB (numeración PDB) a columnas del MSA."""
    target_msa_seq = str(alignment[0].seq)
    msa_non_gap_to_col = []
    for col_idx, aa in enumerate(target_msa_seq):
        if aa != gap_char:
            msa_non_gap_to_col.append(col_idx)
    target_no_gap = "".join(aa for aa in target_msa_seq if aa != gap_char)

    residues, pdb_seq, chain_id_used = extract_chain_sequence(structure, chain_id=chain_id)
    if not pdb_seq or not target_no_gap:
        return {}, residues, chain_id_used

    alignments = pairwise2.align.globalms(pdb_seq, target_no_gap, 1, -1, -10, -0.5)
    if not alignments:
        return {}, residues, chain_id_used
    aligned_pdb, aligned_target, *_ = alignments[0]

    mapping = {}
    pdb_idx = 0
    target_idx = 0
    for pdb_char, target_char in zip(aligned_pdb, aligned_target):
        if pdb_char != gap_char and target_char != gap_char:
            msa_col = msa_non_gap_to_col[target_idx]
            resseq = residues[pdb_idx].id[1]
            if resseq not in mapping:
                mapping[resseq] = msa_col
            pdb_idx += 1
            target_idx += 1
        elif pdb_char != gap_char and target_char == gap_char:
            pdb_idx += 1
        elif pdb_char == gap_char and target_char != gap_char:
            target_idx += 1

    return mapping, residues, chain_id_used
from logging_utils import (
    confirm_file,
    ensure_parent_dir,
    log_info,
    require_file,
)


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


def main():
    pdb_input = snakemake.input.pdb
    msa_input = snakemake.input.msa
    pdb_output = snakemake.output.pdb_conserved
    chain_id = getattr(snakemake.params, "chain", None)

    print("🧮 Calculando conservación evolutiva...")
    alignment = AlignIO.read(msa_input, "fasta")
    scores = calculate_shannon_entropy(alignment)

    # Mapear al PDB usando alineamiento global PDB <-> MSA
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("Target", pdb_input)
    mapping, pdb_residues, chain_id_used = map_pdb_to_msa(
        alignment,
        structure,
        chain_id=chain_id,
    )

    if not mapping:
        print(
            "⚠️ Advertencia: no se pudo generar un mapa PDB->MSA. "
            "Revisa la cadena o la consistencia del MSA."
        )

    print("🎨 Inyectando puntuaciones en el factor B del PDB...")
    unmapped = 0
    for residue in pdb_residues:
        resseq = residue.id[1]
        if resseq not in mapping:
            unmapped += 1
            continue
        score = scores[mapping[resseq]]
        for atom in residue:
            atom.set_bfactor(score)

    if unmapped:
        print(
            "⚠️ Advertencia: "
            f"{unmapped} residuos en la cadena {chain_id_used} "
            "no tienen correspondencia en el MSA (posibles gaps)."
        )

    # Guardar PDB con datos inyectados
    io = PDBIO()
    io.set_structure(structure)
    io.save(pdb_output)

    print(f"✅ ¡Éxito! Archivo generado: {pdb_output}")
    print("   -> Abre este archivo en PyMOL y colorea por B-factor para ver la conservación.")


if __name__ == "__main__":
    main()
require_file(pdb_input, "PDB de entrada")
require_file(msa_input, "MSA de entrada")
ensure_parent_dir(pdb_output)

log_info("🧮 Calculando conservación evolutiva...")
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

log_info("🎨 Inyectando puntuaciones en el factor B del PDB...")

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

confirm_file(pdb_output, "PDB con conservación")
log_info(f"✅ ¡Éxito! Archivo generado: {pdb_output}")
log_info("-> Abre este archivo en PyMOL y colorea por B-factor para ver la conservación.")
