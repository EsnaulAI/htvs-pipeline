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
        input=SimpleNamespace(smi=config["chemistry"]["filtered_library"]),
        output=SimpleNamespace(enumerated="results/module2/library_enumerated.smi"),
        params=SimpleNamespace(),
        wildcards=SimpleNamespace(),
    )

from collections import defaultdict
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdChemReactions

from logging_utils import (
    confirm_file,
    ensure_parent_dir,
    log_info,
    log_warn,
    require_file,
)

PH_MIN = 6.4
PH_MAX = 8.4


def try_dimorphite(smiles):
    try:
        from dimorphite_dl import DimorphiteDL
    except ImportError:
        return None

    dimorphite = DimorphiteDL(
        min_ph=PH_MIN,
        max_ph=PH_MAX,
        label_states=False,
    )

    for attr in ("protonate", "enumerate"):
        if hasattr(dimorphite, attr):
            result = getattr(dimorphite, attr)(smiles)
            return list(result) if result is not None else []

    if callable(dimorphite):
        result = dimorphite(smiles)
        return list(result) if result is not None else []

    raise AttributeError("DimorphiteDL no expone un método compatible para enumeración.")


def canonical_smiles(mol):
    return Chem.MolToSmiles(mol, isomericSmiles=True)


def apply_reaction(mol, reaction):
    products = []
    try:
        for product_tuple in reaction.RunReactants((mol,)):
            product = product_tuple[0]
            Chem.SanitizeMol(product)
            products.append(product)
    except Exception:
        return []
    return products


def enumerate_with_rdkit(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []

    reactions = [
        rdChemReactions.ReactionFromSmarts("[CX3](=O)[OX2H1:1]>>[CX3](=O)[O-:1]"),
        rdChemReactions.ReactionFromSmarts("[CX3](=O)[O-:1]>>[CX3](=O)[OX2H1:1]"),
        rdChemReactions.ReactionFromSmarts("[NX3;H2;!$(NC=O):1]>>[NX4+;H3;!$(NC=O):1]"),
        rdChemReactions.ReactionFromSmarts("[NX3;H1;!$(NC=O):1]>>[NX4+;H2;!$(NC=O):1]"),
        rdChemReactions.ReactionFromSmarts("[NX3;H0;!$(NC=O):1]>>[NX4+;H1;!$(NC=O):1]"),
        rdChemReactions.ReactionFromSmarts("[NX4+;H3,H2,H1;!$(NC=O):1]>>[NX3;H2,H1,H0;!$(NC=O):1]"),
    ]

    seen = {canonical_smiles(mol)}
    queue = [mol]

    for reaction in reactions:
        next_queue = []
        for current in queue:
            for product in apply_reaction(current, reaction):
                smiles_out = canonical_smiles(product)
                if smiles_out not in seen:
                    seen.add(smiles_out)
                    next_queue.append(product)
        queue.extend(next_queue)

    return sorted(seen)


def main():
    input_file = snakemake.input[0]
    output_file = snakemake.output[0]

    require_file(input_file, "SMILES filtrados")
    ensure_parent_dir(output_file)

    log_info("⚗️ Enumerando protómeros (pH 7.4 ± 1).")

    total = 0
    total_states = 0
    missing = 0
    counts = defaultdict(int)

    with open(output_file, "w") as out_f:
        with open(input_file, "r") as in_f:
            for line in in_f:
                if not line.strip():
                    continue
                parts = line.split()
                smiles = parts[0]
                mol_id = parts[1] if len(parts) > 1 else f"Ligand_{total}"

                total += 1

                states = None
                try:
                    states = try_dimorphite(smiles)
                except Exception:
                    states = None

                if not states:
                    states = enumerate_with_rdkit(smiles)

                if not states:
                    missing += 1
                    log_warn(f"No se pudo enumerar: {mol_id}")
                    continue

                for idx, state in enumerate(states, start=1):
                    state_id = f"{mol_id}_protomer{idx}"
                    out_f.write(f"{state}\t{state_id}\n")
                    total_states += 1

                counts[len(states)] += 1

    confirm_file(output_file, "SMILES enumerados")

    log_info(
        "✅ Protomerización completada: "
        f"{total_states} estados generados para {total} ligandos."
    )
    if missing:
        log_warn(f"Ligandos sin protómeros: {missing}")
    if counts:
        summary = ", ".join(f"{k}→{v}" for k, v in sorted(counts.items()))
        log_info(f"Distribución de estados (n_estados→ligandos): {summary}")


if __name__ == "__main__":
    main()
