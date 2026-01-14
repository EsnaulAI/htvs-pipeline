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
