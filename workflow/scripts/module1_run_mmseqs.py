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
        input=SimpleNamespace(pdb="", fasta=config["structure"]["target_fasta"], xml="", original_fasta="", msa=""),
        output=SimpleNamespace(
            pdb="",
            fasta="",
            xml="",
            fasta_msa="",
            fasta_homologs=config["evolution"]["unaligned_homologs"],
            method_log="results/module1/homolog_method.log",
        ),
        params=SimpleNamespace(
            pdb_id=config["structure"]["pdb_id"],
            chain=config["structure"]["chain_id"],
            n_hits=config["evolution"]["n_homologs"],
            e_val=config["evolution"]["e_value"],
        ),
        threads=4,
        wildcards=SimpleNamespace()
    )
# ----------------------------

# workflow/scripts/module1_run_mmseqs.py
import subprocess
import os
import sys
import shutil
from Bio import SeqIO
from pathlib import Path
import yaml

def load_config():
    config_path = Path(__file__).resolve().parents[2] / "config" / "config.yaml"
    with open(config_path, "r") as f:
        return yaml.safe_load(f)
from logging_utils import (
    confirm_file,
    ensure_parent_dir,
    log_error,
    log_info,
    log_warn,
    require_file,
)

def run_mmseqs_search(query_fasta, out_msa_ready, config, threads, method_log_path=None):
    db_dir = config["evolution"]["mmseqs_db_dir"]

    # --- PARÁMETROS AGRESIVOS (Corrección) ---
    # Sensibilidad máxima (7.5 es muy alto, encuentra homólogos remotos)
    sensitivity = 7.5
    # Permitimos hasta 3000 secuencias para tener estadística robusta
    max_seqs = 3000
    # E-value relajado (permitimos primos lejanos)
    e_value = 100.0
    # Identidad mínima bajísima (queremos diversidad, no clones)
    min_id = 0.05
    # -----------------------------------------

    db_target = os.path.join(db_dir, "uniref50_db")
    tmp_dir = config["evolution"]["mmseqs_tmp_dir"]

    require_file(query_fasta, "FASTA de entrada")
    ensure_parent_dir(out_msa_ready)

    log_info(f"🚀 Iniciando Búsqueda Profunda MMseqs2 (Sensibilidad {sensitivity})...")

    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)
    os.makedirs(tmp_dir, exist_ok=True)
    os.makedirs(db_dir, exist_ok=True)

    def run_cmd(cmd):
        log_info(f"CMD: {cmd}")
        try:
            subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError:
            log_error("❌ Error ejecutando MMseqs2.")
            sys.exit(1)

    # 1. Base de Datos
    if not os.path.exists(db_target + ".dbtype"):
        log_info("⬇️ Descargando UniRef50...")
        run_cmd(f"mmseqs databases UniRef50 {db_target} {tmp_dir}")

    # 2. Búsqueda Iterativa
    query_db = os.path.join(tmp_dir, "queryDB")
    result_db = os.path.join(tmp_dir, "resultDB")

    run_cmd(f"mmseqs createdb {query_fasta} {query_db}")

    # Agregamos --start-sens 1 para ir escalando la búsqueda
    # Agregamos -e {e_value} para capturar todo
    cmd_search = (
        f"mmseqs search {query_db} {db_target} {result_db} {tmp_dir} "
        f"-s {sensitivity} --num-iterations 4 --max-seqs {max_seqs} "
        f"--threads {threads} --min-seq-id {min_id} -e {e_value} "
        f"--cov-mode 0 -c 0.3"  # Cobertura requerida bajada al 30%
    )
    run_cmd(cmd_search)

    # 3. Extracción
    fasta_res = os.path.join(tmp_dir, "results.fasta")
    run_cmd(f"mmseqs result2flat {query_db} {db_target} {result_db} {fasta_res} --use-fasta-header")

    # 4. Filtrado Inteligente
    log_info("QC: Filtrando resultados redundantes o basura...")
    original_rec = SeqIO.read(query_fasta, "fasta")
    found_recs = list(SeqIO.parse(fasta_res, "fasta"))

    unique_recs = []
    seen_seqs = set()

    # Siempre incluimos el target original primero
    unique_recs.append(original_rec)
    seen_seqs.add(str(original_rec.seq))

    for rec in found_recs:
        s_str = str(rec.seq)
        # Evitar duplicados exactos
        if s_str in seen_seqs:
            continue
        # Evitar secuencias rotas (muy cortas)
        if len(rec.seq) < len(original_rec.seq) * 0.4:
            continue

        unique_recs.append(rec)
        seen_seqs.add(s_str)

    log_info(f"📊 Secuencias encontradas: {len(found_recs)}")
    log_info(f"✅ Secuencias únicas válidas: {len(unique_recs)}")

    if len(unique_recs) < 50:
        log_warn("⚠️ ALERTA: Aún tenemos pocas secuencias. La conservación será poco confiable.")

    SeqIO.write(unique_recs, out_msa_ready, "fasta")
    confirm_file(out_msa_ready, "FASTA de homologos")

    if method_log_path:
        ensure_parent_dir(method_log_path)
        with open(method_log_path, "w") as log_handle:
            log_handle.write("mmseqs\n")


def main():
    query_fasta = snakemake.input.fasta
    out_msa_ready = snakemake.output.fasta_homologs
    method_log = getattr(snakemake.output, "method_log", None)
    config = load_config()
    threads = snakemake.threads

    run_mmseqs_search(query_fasta, out_msa_ready, config, threads, method_log_path=method_log)


if __name__ == "__main__":
    main()
