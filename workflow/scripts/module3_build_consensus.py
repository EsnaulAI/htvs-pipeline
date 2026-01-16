import os

import pandas as pd


gnina_file = snakemake.input.gnina
top_candidates_file = snakemake.input.top
final_hits_file = getattr(snakemake.input, "final_hits", None)
output_consensus = snakemake.output.consensus
output_report = snakemake.output.report

print("📊 Generando Consenso Final (Gnina + filtro de contactos)...")

try:
    df_gnina = pd.read_csv(gnina_file)
    df_top = pd.read_csv(top_candidates_file)
    df_hits = None
    if final_hits_file and os.path.exists(final_hits_file):
        df_hits = pd.read_csv(final_hits_file)

    if df_hits is not None and not df_hits.empty:
        if "Ligand" in df_hits.columns and "Ligand" in df_gnina.columns:
            df_base = df_gnina[df_gnina["Ligand"].isin(df_hits["Ligand"])]
        else:
            df_base = df_hits
    else:
        df_base = df_gnina

    if "Ligand" in df_base.columns and "ID" in df_top.columns:
        df_final = pd.merge(df_base, df_top, left_on="Ligand", right_on="ID", how="inner")
    elif "ID" in df_base.columns and "ID" in df_top.columns:
        df_final = pd.merge(df_base, df_top, on="ID", how="inner")
    else:
        df_final = df_base

    if "CNNaffinity" in df_final.columns:
        df_final = df_final.sort_values(by="CNNaffinity", ascending=False)

    df_final.to_csv(output_consensus, index=False)
    print(f"✅ Consenso guardado: {output_consensus}")

    top_3 = df_final.head(3)
    with open(output_report, "w", encoding="utf-8") as handle:
        handle.write("=== REPORTE FINAL GHOST JAMMER ===\n")
        handle.write("Status: Ejecucion Exitosa\n\n")
        handle.write("🏆 TOP 3 CANDIDATOS:\n")
        for index, row in top_3.iterrows():
            handle.write(
                f"1. {row.get('Ligand', 'Unknown')} | Score: {row.get('CNNaffinity', 'N/A')}\n"
            )
except Exception as exc:
    print(f"❌ Error en consenso: {exc}")
    os.makedirs(os.path.dirname(output_consensus), exist_ok=True)
    os.makedirs(os.path.dirname(output_report), exist_ok=True)
    with open(output_consensus, "w", encoding="utf-8") as handle:
        handle.write("Error")
    with open(output_report, "w", encoding="utf-8") as handle:
        handle.write(f"Error: {exc}")
