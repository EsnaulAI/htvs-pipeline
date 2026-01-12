# workflow/scripts/generate_chimerax_scene.py

# --- MOCK PARA DESARROLLO (Pylance no se quejará) ---
if "snakemake" not in globals():
    # Esto simula la variable snakemake para que el editor no marque error
    from types import SimpleNamespace
    snakemake = SimpleNamespace(
        input=SimpleNamespace(pdb="results/module1/target_conserved.pdb"),
        output=SimpleNamespace(cxc="results/module1/view_scene.cxc"),
        params=SimpleNamespace(pdb_id="6OCR", chain="B", target_res=513),
        wildcards=SimpleNamespace()
    )
# ----------------------------------------------------

# Inputs desde Snakemake
pdb_file = snakemake.input.pdb
target_res = int(snakemake.params.target_res)
chain_id = snakemake.params.chain
output_cxc = snakemake.output.cxc

print(f"🎨 Generando escena de ChimeraX para el residuo {target_res}...")

# Tamaño de la caja de docking (Debe coincidir con lo que tienes en config.yaml)
BOX_SIZE = 22 

# Contenido del script .cxc
# ChimeraX usa comandos distintos a PyMOL.
# Asumimos que la conservación está guardada en el B-factor

cxc_content = f"""
# --- Script Automático de Ghost-Jammer para ChimeraX ---

# 1. Cargar Estructura
open {pdb_file}

# 2. Estilo Visual (Calidad de Publicación)
lighting soft
graphics silhouettes true color black width 1.5
style surface

# 3. Colorear por Conservación (B-factor)
# Rango 0 (Variable/Azul) a 4.3 (Conservado/Rojo/Maroon)
color byattribute bfactor palette cyan:white:maroon range 0,4.3

# 4. Resaltar el Sitio Alostérico (Target)
# Seleccionamos el residuo objetivo
select /{chain_id}:{target_res}

# Mostramos sus átomos como esferas doradas
show sel atoms
style sel sphere
color sel gold

# Etiqueta flotante
label sel text "Target {target_res}" height 1.5 color black yoffset 2

# 5. DIBUJAR GRIDBOX (Visualización de la Caja de Docking)
# Crea un cubo de malla amarilla centrado exactamente en el residuo seleccionado
shape box center sel size {BOX_SIZE},{BOX_SIZE},{BOX_SIZE} color yellow mesh true name gridbox

# 6. Enfoque final
view orient
zoom sel 0.8
"""

with open(output_cxc, "w") as f:
    f.write(cxc_content)

print(f"✅ Script guardado en: {output_cxc}")
print("   -> Para ver: Abre ChimeraX y ejecuta 'open results/module1/view_scene.cxc'")