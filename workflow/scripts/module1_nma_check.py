# --- MOCK PARA DESARROLLO ---
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
        input=SimpleNamespace(pdb=config["evolution"]["conservation_pdb"]),
        output=SimpleNamespace(plot=config["analysis"]["nma_profile"]),
        params=SimpleNamespace(
            chain=config["structure"]["chain_id"],
            target_res=config["structure"]["target_residue"],
            target_res_name=config["structure"]["target_residue_name"],
        )
    )
# ----------------------------

from pathlib import Path
import sys
import matplotlib.pyplot as plt
import numpy as np

# Importar explícitamente las funciones de ProDy
from prody import parsePDB, ANM, calcSqFlucts, confProDy
import yaml

# Configuración silenciosa
confProDy(verbosity='none')

def load_config():
    config_path = Path(__file__).resolve().parents[2] / "config" / "config.yaml"
    with open(config_path, "r") as f:
        return yaml.safe_load(f)

def run_nma_analysis(pdb_file, chain_id, target_res_num, target_res_name, output_img):
    print(f"🌊 Iniciando Análisis de Modos Normales (ANM) para {pdb_file}...")
    
    # 1. Cargar estructura
    structure = parsePDB(pdb_file)
    calphas = structure.select(f'chain {chain_id} and calpha')
    
    if calphas is None:
        print(f"❌ Error: No se encontraron Carbonos Alpha en la cadena {chain_id}")
        return

    # 2. Construir Modelo Elástico (ANM)
    anm = ANM('AdeB dynamics')
    anm.buildHessian(calphas)
    anm.calcModes(n_modes=20) 
    
    # 3. Calcular Fluctuaciones
    sq_flucts = calcSqFlucts(anm)
    res_nums = calphas.getResnums()
    
    # Buscar nuestro residuo objetivo
    try:
        target_idx = list(res_nums).index(target_res_num)
        target_fluct = sq_flucts[target_idx]
        print(f"\n🔍 ANÁLISIS DE MOVILIDAD PARA {target_res_name} {target_res_num}:")
        print(f"   Fluctuación Cuadrática: {target_fluct:.4f} Å²")
        
        mean_fluct = np.mean(sq_flucts)
        if target_fluct < mean_fluct:
            print("   CONCLUSIÓN: Es una zona RÍGIDA (Posible ancla o bisagra estática).")
        else:
            print("   CONCLUSIÓN: Es una zona FLEXIBLE (Posible loop o puerta).")
            
    except ValueError:
        print(f"⚠️ El residuo {target_res_num} no está en la selección analizada.")

    # 4. Generar Gráfica
    plt.figure(figsize=(10, 5))
    plt.plot(res_nums, sq_flucts, label='Movilidad Teórica (ANM)', color='black', linewidth=1)
    
    # Marcar nuestro target
    plt.axvline(x=target_res_num, color='red', linestyle='--', label=f'{target_res_name} {target_res_num}')
    
    plt.title('Perfil de Flexibilidad de AdeB (Análisis de Modos Normales)')
    plt.xlabel('Número de Residuo')
    plt.ylabel('Fluctuación Cuadrática (Å²)')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Guardar gráfica
    plt.savefig(output_img, dpi=300)
    print(f"\n📊 Gráfica guardada en: {output_img}")

if __name__ == "__main__":
    # Obtener parámetros desde Snakemake
    pdb_in = snakemake.input.pdb
    chain_in = snakemake.params.chain
    
    # --- CAMBIO CRÍTICO: Recibir variable desde Config ---
    res_in = int(snakemake.params.target_res)
    # ---------------------------------------------------
    config = load_config()
    target_res_name = snakemake.params.target_res_name if hasattr(snakemake.params, "target_res_name") else config["structure"]["target_residue_name"]
    out_plot = snakemake.output.plot
    
    run_nma_analysis(pdb_in, chain_in, res_in, target_res_name, out_plot)
