# Dependencias y herramientas externas

## Entornos Conda (versiones declaradas)
| Componente | Versión/Restricción | Propósito |
| --- | --- | --- |
| Python | 3.10 | Lenguaje base del pipeline principal (`environment.yaml`). |
| Snakemake | sin pin (resolver conda) | Orquestación del flujo de trabajo. |
| pandas | sin pin (resolver conda) | Manipulación de datos tabulares. |
| Biopython | sin pin (resolver conda) | Utilidades bioinformáticas (BLAST, secuencias). |
| RDKit | sin pin (resolver conda) | Química computacional y filtrado de ligandos. |
| JupyterLab | sin pin (resolver conda) | Trabajo interactivo en notebooks. |
| Open Babel | sin pin (resolver conda) | Conversión/limpieza de formatos químicos. |

## Entorno específico PyDCA (workflow/envs/pydca.yaml)
| Componente | Versión/Restricción | Propósito |
| --- | --- | --- |
| Python | 3.8 | Compatibilidad con librerías legacy de PyDCA. |
| NumPy | <1.18 | Compatibilidad con SciPy legacy y PyDCA. |
| SciPy | <1.4 | Compatibilidad con PyDCA legacy. |
| Biopython | sin pin (resolver conda) | Utilidades bioinformáticas. |
| Matplotlib | sin pin (resolver conda) | Gráficas de análisis. |
| PyDCA | sin pin (pip) | Análisis de acoplamiento directo (DCA). |

## Herramientas externas (instalación fuera de conda si aplica)
| Herramienta | Versión objetivo | Propósito |
| --- | --- | --- |
| MMseqs2 | >=13.45111 | Búsqueda profunda de homólogos para MSA. |
| fpocket | >=4.0 | Detección de bolsillos en estructuras PDB. |
| UCSF ChimeraX | >=1.6 | Visualización molecular y escenas (`.cxc`). |
| NCBI BLAST (remoto) | N/A (servicio externo) | Búsquedas BLAST remotas vía `NCBIWWW`. |
| ChEMBL API | N/A (servicio externo) | Minería de fármacos en ChEMBL. |

> Nota: Para herramientas sin pin explícito en los YAML, se recomienda registrar la versión efectiva con `tool --version` en el entorno de ejecución.
# Dependencias mínimas

Estas son las versiones mínimas recomendadas para ejecutar el pipeline.

## Entorno base (Python + Snakemake)

- python >= 3.10
- snakemake >= 7.32.4
- pandas >= 1.5.0
- biopython >= 1.79
- rdkit >= 2023.03.1
- jupyterlab >= 4.0.0
- openbabel >= 3.1.1

## Herramientas por módulo

- MAFFT >= 7.505
- fpocket >= 4.0
- Open Babel >= 3.1.1
- MMseqs2 >= 14.7e284
- PyDCA >= 1.0.0

