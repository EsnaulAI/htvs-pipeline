# --- MOCK PARA DESARROLLO (Pylance no se quejará) ---
if "snakemake" not in globals():
    # Esto simula la variable snakemake para que el editor no marque error
    # y tú sepas qué estructura tiene. Solo sirve mientras editas.
    from types import SimpleNamespace
    snakemake = SimpleNamespace(
        input=SimpleNamespace(pdb="", fasta="", xml="", original_fasta="", msa=""),
        output=SimpleNamespace(
            pdb="", fasta="", xml="", fasta_msa="", pdb_conserved="", pdb_msa_map=""
        ),
        params=SimpleNamespace(pdb_id="7KGY", chain="B", n_hits=10, e_val=0.001),
        wildcards=SimpleNamespace()
    )
# ----------------------------------------------------
import json
import math
from Bio import AlignIO
from Bio.Align import PairwiseAligner
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.Polypeptide import is_aa
from Bio.SeqUtils import seq1

def extract_chain_sequence(structure, chain_id=None):
    model = next(structure.get_models())
    if chain_id:
        if chain_id not in model:
            raise ValueError(f"Cadena '{chain_id}' no encontrada en el PDB.")
        chain = model[chain_id]
    else:
        chain = next(model.get_chains())
    residues = [res for res in chain.get_residues() if is_aa(res, standard=True)]
    sequence = "".join(seq1(res.get_resname(), custom_map={"SEC": "U", "PYL": "O"}) for res in residues)
    return residues, sequence, chain.id


def build_gapped_sequences(alignment, gap_char="-"):
    target = alignment.target
    query = alignment.query
    target_aligned = []
    query_aligned = []
    target_pos = 0
    query_pos = 0

    for (t_start, t_end), (q_start, q_end) in zip(
        alignment.aligned[0], alignment.aligned[1]
    ):
        if t_start > target_pos:
            target_aligned.append(target[target_pos:t_start])
            query_aligned.append(gap_char * (t_start - target_pos))
        if q_start > query_pos:
            target_aligned.append(gap_char * (q_start - query_pos))
            query_aligned.append(query[query_pos:q_start])

        target_aligned.append(target[t_start:t_end])
        query_aligned.append(query[q_start:q_end])
        target_pos = t_end
        query_pos = q_end

    if target_pos < len(target):
        target_aligned.append(target[target_pos:])
        query_aligned.append(gap_char * (len(target) - target_pos))
    if query_pos < len(query):
        target_aligned.append(gap_char * (len(query) - query_pos))
        query_aligned.append(query[query_pos:])

    return "".join(target_aligned), "".join(query_aligned)


def map_pdb_to_msa(alignment, structure, chain_id=None, gap_char="-"):
    """Mapea residuos del PDB (numeración PDB) a columnas del MSA (0-based)."""
    target_msa_seq = str(alignment[0].seq)
    msa_non_gap_to_col = []
    for col_idx, aa in enumerate(target_msa_seq):
        if aa != gap_char:
            msa_non_gap_to_col.append(col_idx)
    target_no_gap = "".join(aa for aa in target_msa_seq if aa != gap_char)

    residues, pdb_seq, chain_id_used = extract_chain_sequence(structure, chain_id=chain_id)
    if not pdb_seq or not target_no_gap:
        return {}, residues, chain_id_used

    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    alignments = aligner.align(pdb_seq, target_no_gap)
    alignment = next(iter(alignments), None)
    if alignment is None:
        return {}, residues, chain_id_used
    aligned_pdb, aligned_target = build_gapped_sequences(alignment, gap_char=gap_char)

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


def export_pdb_to_msa_map(mapping, chain_id_used, output_path, msa_length):
    payload = {
        "chain_id": chain_id_used,
        "index_base": 0,
        "msa_length": msa_length,
        "mapping": {str(k): v for k, v in mapping.items()},
    }
    with open(output_path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)

from logging_utils import confirm_file, ensure_parent_dir, log_info, log_warn, require_file


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


pdb_input = snakemake.input.pdb
msa_input = snakemake.input.msa
pdb_output = snakemake.output.pdb_conserved
map_output = getattr(snakemake.output, "pdb_msa_map", None)
chain_id = getattr(snakemake.params, "chain", None)

require_file(pdb_input, "PDB de entrada")
require_file(msa_input, "MSA de entrada")
ensure_parent_dir(pdb_output)
if map_output:
    ensure_parent_dir(map_output)

log_info("🧮 Calculando conservación evolutiva...")
alignment = AlignIO.read(msa_input, "fasta")
scores = calculate_shannon_entropy(alignment)

parser = PDBParser(QUIET=True)
structure = parser.get_structure("Target", pdb_input)
mapping, pdb_residues, chain_id_used = map_pdb_to_msa(
    alignment,
    structure,
    chain_id=chain_id,
)

if not mapping:
    log_warn(
        "⚠️ Advertencia: no se pudo generar un mapa PDB->MSA. "
        "Revisa la cadena o la consistencia del MSA."
    )

if map_output:
    export_pdb_to_msa_map(
        mapping,
        chain_id_used,
        map_output,
        alignment.get_alignment_length(),
    )
    log_info(
        "🧾 Mapa PDB→MSA guardado con índice base 0 (columnas MSA 0-based)."
    )

log_info("🎨 Inyectando puntuaciones en el factor B del PDB...")
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
    log_warn(
        "⚠️ Advertencia: "
        f"{unmapped} residuos en la cadena {chain_id_used} "
        "no tienen correspondencia en el MSA (posibles gaps)."
    )

io = PDBIO()
io.set_structure(structure)
io.save(pdb_output)

confirm_file(pdb_output, "PDB con conservación")
log_info(f"✅ ¡Éxito! Archivo generado: {pdb_output}")
log_info("-> Abre este archivo en PyMOL y colorea por B-factor para ver la conservación.")
