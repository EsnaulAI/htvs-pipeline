# workflow/scripts/module2_filter_entry.py
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import yaml
from rdkit import Chem
from rdkit.Chem import Descriptors

from logging_utils import (
    confirm_file,
    ensure_parent_dir,
    log_info,
    require_file,
)


# --- MOCK PARA DESARROLLO (Pylance no se quejará) ---
if "snakemake" not in globals():
    # Esto simula la variable snakemake para que el editor no marque error
    # y tú sepas qué estructura tiene. Solo sirve mientras editas.

    def load_config():
        config_path = Path(__file__).resolve().parents[2] / "config" / "config.yaml"
        with open(config_path, "r") as f:
            return yaml.safe_load(f)

    config = load_config()
    snakemake = SimpleNamespace(
        input=SimpleNamespace(pdb="", fasta="", xml="", original_fasta="", msa=""),
        output=SimpleNamespace(
            pdb="",
            fasta="",
            xml="",
            fasta_msa="",
            pdb_conserved="",
            filtered_library=config["chemistry"]["filtered_library"],
            deduped_library="results/module2/library_dedup.smi",
        ),
        params=SimpleNamespace(
            pdb_id=config["structure"]["pdb_id"],
            chain=config["structure"]["chain_id"],
            n_hits=config["evolution"]["n_homologs"],
            e_val=config["evolution"]["e_value"],
            ligand_id_prefix=config["chemistry"]["ligand_id_prefix"],
            mw=600,
            logp=5.0,
            amine=True,
        ),
        wildcards=SimpleNamespace(),
    )


def load_config():
    config_path = Path(__file__).resolve().parents[2] / "config" / "config.yaml"
    with open(config_path, "r") as f:
        return yaml.safe_load(f)


# Definir patrón de Amina Primaria (Nitrógeno con 2 hidrógenos, no amida)
PRIMARY_AMINE = Chem.MolFromSmarts("[NX3;H2;!$(NC=O)]")


def passes_entry_rules(mol, mw_max, logp_max, must_have_amine):
    if mol is None:
        return False, "Invalid Molecule"

    # 1. Peso Molecular (< 600 Da)
    mw = Descriptors.MolWt(mol)
    if mw > mw_max:
        return False, f"MW High ({mw:.1f})"

    # 2. LogP
    logp = Descriptors.MolLogP(mol)
    if logp > logp_max:
        return False, f"LogP High ({logp:.2f})"

    # 3. Amina Primaria (Regla clave para Gram-negativas)
    if must_have_amine and not mol.HasSubstructMatch(PRIMARY_AMINE):
        return False, "No Primary Amine"

    return True, "Pass"


def main():
    input_file = snakemake.input[0]
    output_smi = snakemake.output[0]
    output_csv = snakemake.output[1]
    output_dedup = snakemake.output[2]
    mw_max = snakemake.params.mw
    logp_max = snakemake.params.logp
    must_have_amine = snakemake.params.amine

    require_file(input_file, "SMILES de entrada")
    ensure_parent_dir(output_smi)
    ensure_parent_dir(output_csv)
    ensure_parent_dir(output_dedup)

    log_info(f"🧪 Iniciando filtrado eNTRy sobre: {input_file}")

    passed_count = 0
    total_count = 0
    deduped_count = 0
    report_data = []
    seen_inchikeys = set()

    # Abrir archivos
    with open(output_smi, "w") as out_f, open(output_dedup, "w") as dedup_f:
        # Leer SMILES (asumiendo formato: SMILES ID)
        with open(input_file, "r") as in_f:
            for line in in_f:
                if not line.strip():
                    continue
                parts = line.split()
                smiles = parts[0]
                mol_id = parts[1] if len(parts) > 1 else f"Ligand_{total_count}"

                mol = Chem.MolFromSmiles(smiles)
                total_count += 1

                is_valid, reason = passes_entry_rules(mol, mw_max, logp_max, must_have_amine)
                inchikey = Chem.MolToInchiKey(mol) if mol is not None else None

                report_data.append(
                    {
                        "ID": mol_id,
                        "Status": "Accepted" if is_valid else "Rejected",
                        "Reason": reason,
                        "SMILES": smiles,
                        "InChIKey": inchikey or "",
                    }
                )

                if is_valid:
                    passed_count += 1
                    out_f.write(f"{smiles}\t{mol_id}\n")
                    if inchikey and inchikey not in seen_inchikeys:
                        seen_inchikeys.add(inchikey)
                        deduped_count += 1
                        dedup_f.write(f"{smiles}\t{mol_id}\n")

    # Guardar reporte
    df = pd.DataFrame(report_data)
    df.to_csv(output_csv, index=False)

    confirm_file(output_smi, "SMILES filtrados")
    confirm_file(output_dedup, "SMILES deduplicados")
    confirm_file(output_csv, "reporte eNTRy")

    log_info(
        "✅ Filtrado completado: "
        f"{passed_count}/{total_count} ligandos aceptados, "
        f"{deduped_count} únicos por InChIKey."
    )


if __name__ == "__main__":
    main()
