# --- MOCK PARA DESARROLLO ---
if "snakemake" not in globals():
    from types import SimpleNamespace
    snakemake = SimpleNamespace(
        input=SimpleNamespace(pdb="results/module1/target_conserved.pdb"),
        output=SimpleNamespace(plot="results/module1/nma_mobility_profile.png"),
        params=SimpleNamespace(chain="B", target_res=513)
    )
# ----------------------------

import sys
import matplotlib.pyplot as plt
import numpy as np

# Importar explícitamente las funciones de ProDy
from prody import parsePDB, ANM, calcSqFlucts, confProDy
from logging_utils import (
    confirm_file,
    ensure_parent_dir,
    log_error,
    log_info,
    log_warn,
    require_file,
)

# Configuración silenciosa
confProDy(verbosity='none')

def run_nma_analysis(pdb_file, chain_id, target_res_num, output_img):
    log_info(f"🌊 Iniciando Análisis de Modos Normales (ANM) para {pdb_file}...")
    
    # 1. Cargar estructura
    structure = parsePDB(pdb_file)
    calphas = structure.select(f'chain {chain_id} and calpha')
    
    if calphas is None:
        log_error(f"❌ Error: No se encontraron Carbonos Alpha en la cadena {chain_id}")
        sys.exit(1)

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
        log_info(f"🔍 ANÁLISIS DE MOVILIDAD PARA GLU {target_res_num}:")
        log_info(f"Fluctuación Cuadrática: {target_fluct:.4f} Å²")
        
        mean_fluct = np.mean(sq_flucts)
        if target_fluct < mean_fluct:
            log_info("CONCLUSIÓN: Es una zona RÍGIDA (Posible ancla o bisagra estática).")
        else:
            log_info("CONCLUSIÓN: Es una zona FLEXIBLE (Posible loop o puerta).")
            
    except ValueError:
        log_warn(f"⚠️ El residuo {target_res_num} no está en la selección analizada.")

    # 4. Generar Gráfica
    plt.figure(figsize=(10, 5))
    plt.plot(res_nums, sq_flucts, label='Movilidad Teórica (ANM)', color='black', linewidth=1)
    
    # Marcar nuestro target
    plt.axvline(x=target_res_num, color='red', linestyle='--', label=f'GLU {target_res_num}')
    
    plt.title('Perfil de Flexibilidad de AdeB (Análisis de Modos Normales)')
    plt.xlabel('Número de Residuo')
    plt.ylabel('Fluctuación Cuadrática (Å²)')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Guardar gráfica
    plt.savefig(output_img, dpi=300)
    confirm_file(output_img, "gráfica NMA")
    log_info(f"📊 Gráfica guardada en: {output_img}")

if __name__ == "__main__":
    # Obtener parámetros desde Snakemake
    pdb_in = snakemake.input.pdb
    chain_in = snakemake.params.chain
    
    # --- CAMBIO CRÍTICO: Recibir variable desde Config ---
    res_in = int(snakemake.params.target_res)
    # ---------------------------------------------------
    
    out_plot = snakemake.output.plot
    
    require_file(pdb_in, "PDB de entrada")
    ensure_parent_dir(out_plot)

    run_nma_analysis(pdb_in, chain_in, res_in, out_plot)
