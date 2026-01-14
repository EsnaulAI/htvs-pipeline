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
# workflow/scripts/module2_prep_3d.py
from rdkit import Chem
from rdkit.Chem import AllChem
from concurrent.futures import ProcessPoolExecutor
import sys
import os
from logging_utils import (
    confirm_file,
    ensure_parent_dir,
    log_info,
    require_file,
)

def process_ligand(line):
    parts = line.strip().split()
    if not parts: return None
    
    smi = parts[0]
    # Guardamos el nombre en una variable de Python simple para que no se pierda
    name_str = parts[1] if len(parts) > 1 else "Unknown"
    
    try:
        mol = Chem.MolFromSmiles(smi)
        if mol is None: return None
        
        # 1. Hidrógenos a pH 7.4
        mol = Chem.AddHs(mol)
        
        # 2. Generar Múltiples Confórmeros
        params = AllChem.ETKDGv3()
        params.useSmallRingTorsions = True
        # Intentamos generar 5 conformaciones
        ids = AllChem.EmbedMultipleConfs(mol, numConfs=5, params=params)
        
        if not ids:
            ids = AllChem.EmbedMultipleConfs(mol, numConfs=5, useRandomCoords=True)
            if not ids: return None
            
        # 3. Minimización de Energía
        best_energy = float('inf')
        best_conf_id = -1
        
        # Optimizar todas
        results = AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=0)
        
        for i, (converged, energy) in enumerate(results):
            if energy < best_energy:
                best_energy = energy
                best_conf_id = i
                
        if best_conf_id == -1: return None
        
        # 4. Crear la molécula final limpia
        final_mol = Chem.Mol(mol)
        final_mol.RemoveAllConformers()
        final_mol.AddConformer(mol.GetConformer(best_conf_id))
        
        # --- CORRECCIÓN CRÍTICA: Re-asignar el nombre explícitamente ---
        # Aseguramos que la propiedad _Name viaje con el objeto final
        final_mol.SetProp("_Name", name_str)
        
        return final_mol
    except Exception:
        return None

def main():
    input_smi = snakemake.input.smi
    out_dir = snakemake.params.out_dir
    threads = snakemake.threads
    output_flag = snakemake.output[0] # El archivo .flag

    require_file(input_smi, "SMILES de entrada")
    ensure_parent_dir(output_flag)

    log_info(f"⚡ TURBO V2 (FIXED): Procesando ligandos con {threads} núcleos...")
    os.makedirs(out_dir, exist_ok=True)
    
    with open(input_smi, 'r') as f:
        lines = f.readlines()
    
    count = 0
    
    # Paralelización
    with ProcessPoolExecutor(max_workers=threads) as executor:
        for mol in executor.map(process_ligand, lines):
            if mol:
                # --- BLINDAJE EXTRA: Verificar si trae nombre ---
                if mol.HasProp("_Name"):
                    name = mol.GetProp("_Name")
                else:
                    # Si por alguna razón se perdió, inventamos uno para no romper el pipeline
                    name = f"Ligand_{count}_recovered"
                
                # Limpiar nombre de caracteres peligrosos para archivos
                safe_name = "".join([c for c in name if c.isalnum() or c in ('_','-')])
                if not safe_name: safe_name = f"Ligand_{count}"
                
                outfile = os.path.join(out_dir, f"{safe_name}.sdf")
                
                w = Chem.SDWriter(outfile)
                w.write(mol)
                w.close()
                count += 1

    # Crear el flag para que Snakemake sepa que terminamos
    with open(output_flag, 'w') as f:
        f.write("Done")

    confirm_file(output_flag, "flag de salida")
    log_info(f"✅ Procesamiento Terminado: {count} ligandos listos en {out_dir}")

if __name__ == "__main__":
    main()
