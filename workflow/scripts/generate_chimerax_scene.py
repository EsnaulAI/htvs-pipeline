# workflow/scripts/generate_chimerax_scene.py

# --- MOCK PARA DESARROLLO (Pylance no se quejará) ---
if "snakemake" not in globals():
    # Esto simula la variable snakemake para que el editor no marque error
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
        output=SimpleNamespace(cxc=config["visualization"]["view_scene"]),
        params=SimpleNamespace(
            pdb_id=config["structure"]["pdb_id"],
            chain=config["structure"]["chain_id"],
            target_res=config["structure"]["target_residue"],
            target_res_name=config["structure"]["target_residue_name"],
        ),
        wildcards=SimpleNamespace()
    )
# ----------------------------------------------------

from pathlib import Path
import yaml

def load_config():
    config_path = Path(__file__).resolve().parents[2] / "config" / "config.yaml"
    with open(config_path, "r") as f:
        return yaml.safe_load(f)

# Inputs desde Snakemake
pdb_file = snakemake.input.pdb
target_res = int(snakemake.params.target_res)
chain_id = snakemake.params.chain
output_cxc = snakemake.output.cxc
config = load_config()
target_res_name = snakemake.params.target_res_name if hasattr(snakemake.params, "target_res_name") else config["structure"]["target_residue_name"]

print(f"🎨 Generando escena de ChimeraX para el residuo {target_res_name} {target_res}...")

# Tamaño de la caja de docking (Debe coincidir con lo que tienes en config.yaml)
BOX_SIZE = config["docking"]["size_x"]

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
label sel text "Target {target_res_name} {target_res}" height 1.5 color black yoffset 2

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
print(f"   -> Para ver: Abre ChimeraX y ejecuta 'open {output_cxc}'")
