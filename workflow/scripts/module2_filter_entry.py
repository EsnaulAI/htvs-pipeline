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
# workflow/scripts/module2_filter_entry.py
from rdkit import Chem
from rdkit.Chem import Descriptors
import pandas as pd
import sys

# Definir patrón de Amina Primaria (Nitrógeno con 2 hidrógenos, no amida)
PRIMARY_AMINE = Chem.MolFromSmarts("[NX3;H2;!$(NC=O)]")

def passes_entry_rules(mol):
    if mol is None: return False, "Invalid Molecule"
    
    # 1. Peso Molecular (< 600 Da)
    mw = Descriptors.MolWt(mol)
    if mw > 600:
        return False, f"MW High ({mw:.1f})"
        
    # 2. Amina Primaria (Regla clave para Gram-negativas)
    if not mol.HasSubstructMatch(PRIMARY_AMINE):
        return False, "No Primary Amine"
        
    return True, "Pass"

def main():
    input_file = snakemake.input[0]
    output_smi = snakemake.output[0]
    output_csv = snakemake.output[1]
    
    print(f"🧪 Iniciando filtrado eNTRy sobre: {input_file}")
    
    passed_count = 0
    total_count = 0
    report_data = []
    
    # Abrir archivos
    with open(output_smi, 'w') as out_f:
        # Leer SMILES (asumiendo formato: SMILES ID)
        with open(input_file, 'r') as in_f:
            for line in in_f:
                if not line.strip(): continue
                parts = line.split()
                smiles = parts[0]
                mol_id = parts[1] if len(parts) > 1 else f"Ligand_{total_count}"
                
                mol = Chem.MolFromSmiles(smiles)
                total_count += 1
                
                is_valid, reason = passes_entry_rules(mol)
                
                report_data.append({
                    "ID": mol_id,
                    "Status": "Accepted" if is_valid else "Rejected",
                    "Reason": reason,
                    "SMILES": smiles
                })
                
                if is_valid:
                    passed_count += 1
                    out_f.write(f"{smiles}\t{mol_id}\n")

    # Guardar reporte
    df = pd.DataFrame(report_data)
    df.to_csv(output_csv, index=False)
    
    print(f"✅ Filtrado completado: {passed_count}/{total_count} ligandos aceptados.")

if __name__ == "__main__":
    main()# workflow/scripts/module2_filter_entry.py
from rdkit import Chem
from rdkit.Chem import Descriptors
import pandas as pd
import sys

# Definir patrón de Amina Primaria (Nitrógeno con 2 hidrógenos, no amida)
PRIMARY_AMINE = Chem.MolFromSmarts("[NX3;H2;!$(NC=O)]")

def passes_entry_rules(mol):
    if mol is None: return False, "Invalid Molecule"
    
    # 1. Peso Molecular (< 600 Da)
    mw = Descriptors.MolWt(mol)
    if mw > 600:
        return False, f"MW High ({mw:.1f})"
        
    # 2. Amina Primaria (Regla clave para Gram-negativas)
    if not mol.HasSubstructMatch(PRIMARY_AMINE):
        return False, "No Primary Amine"
        
    return True, "Pass"

def main():
    input_file = snakemake.input[0]
    output_smi = snakemake.output[0]
    output_csv = snakemake.output[1]
    
    print(f"🧪 Iniciando filtrado eNTRy sobre: {input_file}")
    
    passed_count = 0
    total_count = 0
    report_data = []
    
    # Abrir archivos
    with open(output_smi, 'w') as out_f:
        # Leer SMILES (asumiendo formato: SMILES ID)
        with open(input_file, 'r') as in_f:
            for line in in_f:
                if not line.strip(): continue
                parts = line.split()
                smiles = parts[0]
                mol_id = parts[1] if len(parts) > 1 else f"Ligand_{total_count}"
                
                mol = Chem.MolFromSmiles(smiles)
                total_count += 1
                
                is_valid, reason = passes_entry_rules(mol)
                
                report_data.append({
                    "ID": mol_id,
                    "Status": "Accepted" if is_valid else "Rejected",
                    "Reason": reason,
                    "SMILES": smiles
                })
                
                if is_valid:
                    passed_count += 1
                    out_f.write(f"{smiles}\t{mol_id}\n")

    # Guardar reporte
    df = pd.DataFrame(report_data)
    df.to_csv(output_csv, index=False)
    
    print(f"✅ Filtrado completado: {passed_count}/{total_count} ligandos aceptados.")

if __name__ == "__main__":
    main()