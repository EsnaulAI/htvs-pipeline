# Convenciones de nombres

## Prefijos de módulos

Para asegurar consistencia en el pipeline, todos los archivos, scripts y artefactos derivados deben usar uno de estos prefijos de módulo:

- `module1_`
- `module2_`
- `module3_`

## Rutas de resultados

Los resultados deben almacenarse en rutas estandarizadas por módulo:

- `results/module1/`
- `results/module2/`
- `results/module3/`

## Reglas rápidas

- Los identificadores generados (p. ej., IDs de ligandos o marcas de análisis) deben respetar el prefijo del módulo correspondiente.
- No usar rutas absolutas en scripts; leer rutas desde `config/config.yaml`.
