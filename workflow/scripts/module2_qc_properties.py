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
        output=SimpleNamespace(report="results/module2/qc_report.csv"),
        params=SimpleNamespace(),
        wildcards=SimpleNamespace(),
    )

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors
import numpy as np
import pandas as pd

from logging_utils import (
    confirm_file,
    ensure_parent_dir,
    log_info,
    log_warn,
    require_file,
)


def compute_properties(mol):
    return {
        "MW": Descriptors.MolWt(mol),
        "LogP": Crippen.MolLogP(mol),
        "TPSA": rdMolDescriptors.CalcTPSA(mol),
        "HBD": Lipinski.NumHDonors(mol),
        "HBA": Lipinski.NumHAcceptors(mol),
    }


def main():
    input_file = snakemake.input[0]
    output_csv = snakemake.output[0]

    require_file(input_file, "SMILES filtrados")
    ensure_parent_dir(output_csv)

    log_info(f"📊 Calculando propiedades fisicoquímicas desde: {input_file}")

    properties = {"MW": [], "LogP": [], "TPSA": [], "HBD": [], "HBA": []}
    invalid_count = 0
    total_count = 0

    with open(input_file, "r") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            total_count += 1
            smiles = line.split()[0]
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                invalid_count += 1
                continue
            props = compute_properties(mol)
            for key, value in props.items():
                properties[key].append(value)

    if invalid_count:
        log_warn(f"⚠️ Se omitieron {invalid_count} moléculas inválidas de {total_count}.")

    if not properties["MW"]:
        raise ValueError("No se pudieron calcular propiedades: lista de entrada vacía o inválida.")

    rows = []
    for key, values in properties.items():
        p5, p50, p95 = np.percentile(values, [5, 50, 95])
        rows.append({
            "property": key,
            "p5": round(float(p5), 4),
            "p50": round(float(p50), 4),
            "p95": round(float(p95), 4),
        })

    df = pd.DataFrame(rows)
    df.to_csv(output_csv, index=False)

    confirm_file(output_csv, "reporte de QC")
    log_info(f"✅ Reporte QC generado: {output_csv}")


if __name__ == "__main__":
    main()
