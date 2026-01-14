# Checklist de ejecución

1. **Activar entorno Conda/Mamba**
   - Comando:
     ```bash
     mamba activate htvs-pipeline
     ```
   - Salida esperada (ejemplo):
     ```text
     (htvs-pipeline) $
     ```

2. **Validar `config/config.yaml`**
   - Comando:
     ```bash
     yamllint config/config.yaml
     ```
   - Salida esperada (ejemplo):
     ```text
     config/config.yaml
       0:0       warning  missing document start "---"  (document-start)
     ```

3. **Ejecutar Snakemake con reglas clave**
   - Comando:
     ```bash
     snakemake --cores 4 --use-conda all qc align
     ```
   - Salida esperada (ejemplo):
     ```text
     Building DAG of jobs...
     [Tue Mar 19 10:22:11 2024]
     rule qc:
         input: data/reads/sample_R1.fastq.gz, data/reads/sample_R2.fastq.gz
         output: results/qc/sample_fastqc.html
     ```

4. **Verificar outputs en `results/`**
   - Comando:
     ```bash
     ls -lah results/
     ```
   - Salida esperada (ejemplo):
     ```text
     drwxr-xr-x  6 user user 4.0K Mar 19 10:30 .
     drwxr-xr-x 12 user user 4.0K Mar 19 10:00 ..
     -rw-r--r--  1 user user 1.2M Mar 19 10:30 summary.tsv
     ```
