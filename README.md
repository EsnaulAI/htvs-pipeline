# htvs-pipeline

Pipeline Snakemake para screening virtual y análisis bioinformático.

## Requisitos mínimos

- Python >= 3.10
- Snakemake >= 7.32
- MAFFT >= 7.520
- MMseqs2 >= 15.6f452
- fpocket >= 4.0
- Open Babel >= 3.1.1

## Instalación reproducible (Conda/Mamba)

```bash
mamba env create -f environment.yaml
mamba activate htvs_main
```

Para ejecutar reglas con ambientes por regla:

```bash
snakemake --use-conda --cores 4
```

Si prefieres instalar manualmente herramientas externas, usa los entornos en
`workflow/envs/*.yaml` como plantillas reproducibles.

## Solución de problemas (Conda)

Si Snakemake falla con `Non-conda folder exists at prefix`, limpia el
directorio de cache de entornos y vuelve a intentar:

```bash
rm -rf .snakemake/conda
snakemake --use-conda --cores 4 -s workflow/Snakefile
```

También puedes usar el script auxiliar:

```bash
./scripts/clean_snakemake_conda.sh
snakemake --use-conda --cores 4 -s workflow/Snakefile
```

Alternativamente, puedes usar un prefijo nuevo para los entornos:

```bash
snakemake --use-conda --conda-prefix .snakemake/conda-clean --cores 4 -s workflow/Snakefile
```

Si el error persiste con mamba, prueba el frontend clásico de conda con un
prefijo dedicado (evita directorios corruptos dentro del repo):

```bash
snakemake --use-conda --conda-frontend conda --conda-prefix ~/.snakemake/conda-htvs --cores 4 -s workflow/Snakefile
```

También puedes usar el perfil incluido:

```bash
./scripts/run_snakemake_conda.sh --cores 4 -s workflow/Snakefile
```

Si aparece un error tipo `AttributeError: Can't get attribute 'Namedlist._used_attribute'`,
asegúrate de recrear los entornos con la versión de Snakemake fijada en los
YAMLs (borra el prefijo y relanza):

```bash
rm -rf ~/.snakemake/conda-htvs
./scripts/run_snakemake_conda.sh --cores 4 -s workflow/Snakefile
```

Si aparece un `LockException` indicando que el directorio no puede ser
bloqueado, libera el bloqueo y reintenta:

```bash
snakemake --unlock -s workflow/Snakefile
./scripts/run_snakemake_conda.sh --cores 4 -s workflow/Snakefile
```

También puedes usar el script auxiliar:

```bash
./scripts/unlock_snakemake.sh
./scripts/run_snakemake_conda.sh --cores 4 -s workflow/Snakefile
```

## Gnina vía contenedor (Docker/Singularity)

Si no puedes instalar Gnina con conda/mamba, la regla `run_gnina_rescore`
usa un wrapper que ejecuta Gnina en contenedor. Por defecto intenta Docker
y luego Apptainer/Singularity.

Ejemplo directo (Docker):

```bash
docker run --rm -v "$PWD:/work" -w /work gnina/gnina gnina --help
```

Ejemplo directo (Singularity/Apptainer con imagen local `.sif`):

```bash
apptainer exec --bind "$PWD:/work" /ruta/gnina.sif gnina --help
```

Puedes personalizar la imagen y el runtime usando variables de entorno:

```bash
export GNINA_CONTAINER_IMAGE=gnina/gnina
export GNINA_CONTAINER_RUNTIME=docker  # o apptainer/singularity
```
