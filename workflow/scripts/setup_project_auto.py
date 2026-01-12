import requests
import pandas as pd
import yaml
import sys
import os

# CONFIGURACIÓN
ORGANISM = "Acinetobacter baumannii"
KEYWORDS = ["AdeB", "Efflux", "Pump", "Transporter"] 
CONFIG_FILE = "config/config.yaml"

print(f"🚀 BUSCANDO {ORGANISM} (MODO SIN CENSURA)...")

# 1. BÚSQUEDA AMPLIA
query = {
    "query": {
        "type": "terminal",
        "service": "text",
        "parameters": {
            "attribute": "rcsb_entity_source_organism.taxonomy_lineage.name",
            "operator": "contains_phrase",
            "value": ORGANISM
        }
    },
    "return_type": "entry",
    "request_options": {"return_all_hits": True}
}

try:
    resp = requests.post("https://search.rcsb.org/rcsbsearch/v2/query", json=query)
    pdb_ids = [item["identifier"] for item in resp.json().get("result_set", [])]
    
    print(f"   -> {len(pdb_ids)} estructuras totales encontradas.")

    # 2. DESCARGA Y FILTRADO
    gql_query = """
    query structure_info($ids: [String!]!) {
      entries(entry_ids: $ids) {
        rcsb_id
        struct { title }
        rcsb_entry_info {
          resolution_combined
          structure_determination_methodology
          deposition_date
        }
      }
    }
    """
    
    candidates = []
    chunk_size = 50
    data_url = "https://data.rcsb.org/graphql"

    print("   -> Analizando calidad...")
    for i in range(0, len(pdb_ids), chunk_size):
        chunk = pdb_ids[i:i+chunk_size]
        try:
            r = requests.post(data_url, json={"query": gql_query, "variables": {"ids": chunk}})
            entries = r.json()["data"]["entries"]
            
            for entry in entries:
                if not entry: continue
                
                title = entry["struct"]["title"] or ""
                
                # FILTRO: ¿Suena a bomba de eflujo?
                is_relevant = False
                for k in KEYWORDS:
                    if k.upper() in title.upper():
                        is_relevant = True
                        break
                
                if not is_relevant: continue
                
                # Resolución (Si es nula, ponemos 99.9)
                res_list = entry["rcsb_entry_info"]["resolution_combined"]
                res = res_list[0] if res_list else 99.9
                
                candidates.append({
                    "id": entry["rcsb_id"],
                    "res": res,
                    "title": title[:60],
                    "method": entry["rcsb_entry_info"]["structure_determination_methodology"]
                })
        except: continue

    if not candidates:
        print("❌ NO SE ENCONTRARON BOMBAS AUTOMÁTICAMENTE.")
        # FALLBACK MANUAL: Si todo falla, forzamos la mejor conocida
        winner = {"id": "6OCR", "res": 3.2, "title": "Crystal structure of AdeB (Forzado por sistema)"}
        print(f"⚠️  ACTIVANDO PROTOCOLO DE EMERGENCIA: Usando {winner['id']} (Estándar de Oro).")
    else:
        # Ordenar y ganar
        df = pd.DataFrame(candidates)
        df = df.sort_values(by="res", ascending=True)
        winner = df.iloc[0]
        
        print("\n🏆 TOP 3 CANDIDATOS ENCONTRADOS:")
        print(df.head(3).to_string(index=False))

    # 3. ACTUALIZAR CONFIG
    print(f"\n⚙️  Configurando proyecto con: {winner['id']} ({winner['res']} Å)")
    
    with open(CONFIG_FILE, 'r') as f:
        config = yaml.safe_load(f)
    
    config['structure']['pdb_id'] = winner['id']
    
    with open(CONFIG_FILE, 'w') as f:
        yaml.dump(config, f, sort_keys=False)
        
    print("✅ LISTO. Configuración guardada.")

except Exception as e:
    print(f"Error: {e}")