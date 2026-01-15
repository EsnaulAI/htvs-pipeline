# Estructura del repositorio

## Árbol real (profundidad 2)
```
.
├── README.md
├── config
│   └── config.yaml
├── data
├── docs
├── envs
├── environment.yaml
├── logs
├── models
├── notebooks
├── results
├── scripts
└── workflow
    ├── Snakefile
    ├── envs
    └── scripts
```

## Carpetas estándar
- `docs/`: documentación técnica y descriptiva del proyecto.
- `data/`: datos de entrada y datasets de referencia.
- `results/`: resultados finales del pipeline (tablas, reportes, figuras).
- `workflow/`: definición del flujo de trabajo de Snakemake (Snakefile, envs, scripts específicos).
- `scripts/`: scripts auxiliares y utilitarios fuera del flujo principal.
- `envs/`: archivos de entorno y dependencias reutilizables.
- `logs/`: registros de ejecución y trazas.
- `notebooks/`: cuadernos de análisis exploratorio.

## Extras presentes en el repo
- `config/`: configuración central del pipeline.
- `models/`: modelos entrenados o recursos derivados.

## Rutas críticas (referencias clave)
- `config/config.yaml`: parámetros del proyecto y rutas de trabajo.
- `workflow/Snakefile`: orquestación principal de Snakemake.
- `workflow/scripts/`: scripts ejecutados por reglas del pipeline.
- `workflow/envs/pydca.yaml`: entorno específico para PyDCA.
- `environment.yaml`: entorno general de la ejecución.
- `data/`: insumos y bases de datos locales.
- `results/`: salida del pipeline (reportes, figuras y resultados intermedios/finales).

## Decisión sobre fuente química
- El pipeline descarga la librería base desde ChEMBL usando `workflow/scripts/fetch_chembl_drugs.py`.
- Se eliminó `chemistry.zinc_url` de la configuración porque no se utiliza en el flujo actual; la ruta `chemistry.raw_library` permanece como destino de la librería generada desde ChEMBL.
