"""
Script: generar_tabla_diversidad_alfa.py
Genera las tablas suplementarias S1–S6 (valores individuales por muestra)
y la tabla resumen para el manuscrito (media ± DE por ambiente y comunidad).

Fuentes de datos:
  - 16S: final.opti_mcc.groups.ave-std.summary  (con submuestreo)
  - 18S: final.opti_mcc.groups.summary           (sin submuestreo; conteos bajos)
"""

import pandas as pd
import numpy as np
from pathlib import Path

# ───────────────────────────────────────────────
# RUTAS
# ───────────────────────────────────────────────
BASE_16S = Path("example_data/16S")
BASE_18S = Path("example_data/18S")
OUT_DIR  = Path("example_data/16S")

FILE_16S = BASE_16S / "final.opti_mcc.groups.ave-std.summary"
FILE_18S = BASE_18S / "final.opti_mcc.groups.summary"

# Metadatos: qué tipo es cada muestra
GRUPOS = {
    "AguaP1": "Agua", "AguaP2": "Agua", "AguaP3": "Agua",
    "AguaP4": "Agua", "AguaP5": "Agua",
    "CrassostreavirginicaG1":  "Ostión", "CrassostreavirginicaG2":  "Ostión",
    "CrassostreavirginicaG3":  "Ostión", "CrassostreavirginicaG4":  "Ostión",
    "CrassostreavirginicaG5":  "Ostión", "CrassostreavirginicaG6":  "Ostión",
    "CrassostreavirginicaG7":  "Ostión", "CrassostreavirginicaG8":  "Ostión",
    "CrassostreavirginicaG9":  "Ostión", "CrassostreavirginicaG10": "Ostión",
    "HecesP6":  "Heces humanas", "HecesP7":  "Heces humanas",
    "HecesP8":  "Heces humanas", "HecesP9":  "Heces humanas",
    "HecesP10": "Heces humanas", "HecesP11": "Heces humanas",
    "HecesP12": "Heces humanas", "HecesP13": "Heces humanas",
    "HecesP14": "Heces humanas",
}

# Etiquetas cortas para las muestras en las tablas
ETIQUETAS = {
    "AguaP1": "W1", "AguaP2": "W2", "AguaP3": "W3",
    "AguaP4": "W4", "AguaP5": "W5",
    "CrassostreavirginicaG1":  "O1",  "CrassostreavirginicaG2":  "O2",
    "CrassostreavirginicaG3":  "O3",  "CrassostreavirginicaG4":  "O4",
    "CrassostreavirginicaG5":  "O5",  "CrassostreavirginicaG6":  "O6",
    "CrassostreavirginicaG7":  "O7",  "CrassostreavirginicaG8":  "O8",
    "CrassostreavirginicaG9":  "O9",  "CrassostreavirginicaG10": "O10",
    "HecesP6":  "F1",  "HecesP7":  "F2",  "HecesP8":  "F3",
    "HecesP9":  "F4",  "HecesP10": "F5",  "HecesP11": "F6",
    "HecesP12": "F7",  "HecesP13": "F8",  "HecesP14": "F9",
}

ORDEN_AMB = ["Agua", "Ostión", "Heces humanas"]

# ───────────────────────────────────────────────
# CARGA DE DATOS
# ───────────────────────────────────────────────

def cargar_16S(path: Path) -> pd.DataFrame:
    """Lee el archivo ave-std.summary de mothur y devuelve sólo filas 'ave'."""
    df = pd.read_csv(path, sep="\t")
    df = df[df["method"] == "ave"].copy()
    df = df.rename(columns={"group": "muestra",
                             "chao": "Chao1",
                             "shannon": "Shannon",
                             "invsimpson": "InvSimpson"})
    df["ambiente"] = df["muestra"].map(GRUPOS)
    df["etiqueta"] = df["muestra"].map(ETIQUETAS)
    return df[["muestra", "etiqueta", "ambiente", "Chao1", "Shannon", "InvSimpson"]]


def cargar_18S(path: Path, min_reads: int = 10) -> pd.DataFrame:
    """
    Lee el archivo groups.summary de mothur (sin submuestreo).
    Excluye muestras con valores negativos o con conteos < min_reads.
    El umbral min_reads se aplica según la tabla de conteos conocida.
    """
    # Conteos reales de lecturas en el shared del 18S
    conteos_18S = {
        "AguaP1": 7919, "AguaP2": 680, "AguaP3": 3313,
        "AguaP4": 132, "AguaP5": 742,
        "CrassostreavirginicaG2": 1, "CrassostreavirginicaG3": 2,
        "CrassostreavirginicaG4": 15, "CrassostreavirginicaG5": 1,
        "CrassostreavirginicaG6": 5, "CrassostreavirginicaG7": 4,
        "CrassostreavirginicaG8": 2,
        "HecesP10": 14, "HecesP11": 5, "HecesP12": 39,
        "HecesP13": 302, "HecesP14": 20,
        "HecesP6": 257, "HecesP7": 1, "HecesP8": 5,
        "HecesP9": 2,
    }

    df = pd.read_csv(path, sep="\t")
    df = df.rename(columns={"group": "muestra",
                             "chao": "Chao1",
                             "shannon": "Shannon",
                             "invsimpson": "InvSimpson"})

    # Marcar muestras con datos insuficientes
    df["reads"] = df["muestra"].map(conteos_18S).fillna(0).astype(int)
    df["insuficiente"] = (
        (df["reads"] < min_reads) |
        (df["Shannon"] < 0) |
        (df["InvSimpson"] < 1)
    )

    df["ambiente"] = df["muestra"].map(GRUPOS)
    df["etiqueta"] = df["muestra"].map(ETIQUETAS)
    return df[["muestra", "etiqueta", "ambiente", "Chao1", "Shannon", "InvSimpson",
               "reads", "insuficiente"]]


# ───────────────────────────────────────────────
# TABLAS SUPLEMENTARIAS S1–S6
# ───────────────────────────────────────────────

def tabla_suplementaria(df: pd.DataFrame, ambiente: str, comunidad: str,
                        es_18S: bool = False) -> pd.DataFrame:
    """Genera la tabla individual con valores por muestra para un ambiente dado."""
    sub = df[df["ambiente"] == ambiente].copy()

    if es_18S:
        sub["nota"] = sub["insuficiente"].map(
            {True: "datos insuficientes (< 10 lecturas)", False: ""}
        )
        # Marcar valores no confiables como N/A
        for col in ["Chao1", "Shannon", "InvSimpson"]:
            sub.loc[sub["insuficiente"], col] = np.nan

    cols_out = ["etiqueta", "Chao1", "Shannon", "InvSimpson"]
    if es_18S:
        cols_out.append("nota")

    tabla = sub[cols_out].rename(columns={"etiqueta": "Muestra"})
    tabla = tabla.sort_values("Muestra").reset_index(drop=True)

    # Añadir fila de media y DE al final
    num_cols = ["Chao1", "Shannon", "InvSimpson"]
    media = tabla[num_cols].mean()
    de    = tabla[num_cols].std(ddof=1)
    n_val = tabla[num_cols].notna().sum()

    fila_media = pd.DataFrame([{
        "Muestra": f"Media (n={int(n_val.min())})",
        **{c: round(media[c], 2) for c in num_cols}
    }])
    fila_de = pd.DataFrame([{
        "Muestra": "DE",
        **{c: round(de[c], 2) for c in num_cols}
    }])

    tabla[num_cols] = tabla[num_cols].round(2)
    tabla = pd.concat([tabla, fila_media, fila_de], ignore_index=True)
    return tabla


# ───────────────────────────────────────────────
# TABLAS RESUMEN PARA EL MANUSCRITO
# ───────────────────────────────────────────────

def resumen_por_comunidad(df: pd.DataFrame, comunidad: str,
                          es_18S: bool = False) -> pd.DataFrame:
    """
    Tabla resumen para UNA comunidad (16S o 18S).
    Filas = ambientes, columnas = Chao1, Shannon, InvSimpson (media ± DE).
    """
    filas = []
    for amb in ORDEN_AMB:
        if es_18S:
            sub = df[(df["ambiente"] == amb) & (~df["insuficiente"])]
        else:
            sub = df[df["ambiente"] == amb]

        fila = {"Ambiente": amb}
        for idx in ["Chao1", "Shannon", "InvSimpson"]:
            vals = sub[idx].dropna()
            if len(vals) >= 2:
                fila[idx] = f"{vals.mean():.2f} ± {vals.std(ddof=1):.2f}"
                fila[f"{idx}_n"] = len(vals)
            elif len(vals) == 1:
                fila[idx] = f"{vals.iloc[0]:.2f} (n=1)"
                fila[f"{idx}_n"] = 1
            else:
                fila[idx] = "N/A"
                fila[f"{idx}_n"] = 0
        filas.append(fila)

    df_out = pd.DataFrame(filas)[["Ambiente", "Chao1", "Shannon", "InvSimpson"]]
    return df_out


# ───────────────────────────────────────────────
# MAIN
# ───────────────────────────────────────────────

def main():
    print("Cargando datos 16S...")
    df16 = cargar_16S(FILE_16S)
    print(f"  {len(df16)} muestras cargadas")

    print("Cargando datos 18S...")
    df18 = cargar_18S(FILE_18S, min_reads=10)
    print(f"  {len(df18)} muestras; "
          f"{df18['insuficiente'].sum()} con datos insuficientes (< 10 lecturas)")

    # ── Tablas suplementarias ──────────────────
    config_tablas = [
        ("S1", "Agua",         "16S", df16, False),
        ("S2", "Ostión",       "16S", df16, False),
        ("S3", "Heces humanas","16S", df16, False),
        ("S4", "Agua",         "18S", df18, True),
        ("S5", "Ostión",       "18S", df18, True),
        ("S6", "Heces humanas","18S", df18, True),
    ]

    print("\n=== TABLAS SUPLEMENTARIAS (valores individuales) ===\n")
    tablas_sup = {}
    for etiq, amb, com, df, es18 in config_tablas:
        t = tabla_suplementaria(df, amb, com, es_18S=es18)
        tablas_sup[etiq] = t
        fname = OUT_DIR / f"Tabla_{etiq}_alpha_{com}_{amb.replace(' ', '_')}.csv"
        t.to_csv(fname, index=False, encoding="utf-8-sig")
        print(f"Tabla {etiq} – {com} | {amb}:")
        print(t.to_string(index=False))
        print(f"  → Guardada en {fname}\n")

    # ── Tabla resumen ─────────────────────────
    # ── Tabla resumen combinada 16S + 18S ────────
    print("\n=== TABLA RESUMEN (16S + 18S) ===\n")
    t16 = resumen_por_comunidad(df16, "16S", es_18S=False).add_prefix("16S_").rename(columns={"16S_Ambiente": "Ambiente"})
    t18 = resumen_por_comunidad(df18, "18S", es_18S=True).add_prefix("18S_").rename(columns={"18S_Ambiente": "Ambiente"})
    t_combined = t16.merge(t18, on="Ambiente")
    # Reordenar columnas: Ambiente | 16S Chao1 | 16S Shannon | 16S InvSimpson | 18S Chao1 | 18S Shannon | 18S InvSimpson
    t_combined = t_combined[["Ambiente",
                              "16S_Chao1", "16S_Shannon", "16S_InvSimpson",
                              "18S_Chao1", "18S_Shannon", "18S_InvSimpson"]]
    t_combined.columns = ["Ambiente",
                          "16S Chao1", "16S Shannon", "16S InvSimpson",
                          "18S Chao1", "18S Shannon", "18S InvSimpson"]
    print(t_combined.to_string(index=False))
    f_comb = OUT_DIR / "Tabla_resumen_alpha_diversidad.csv"
    t_combined.to_csv(f_comb, index=False, encoding="utf-8-sig")
    print(f"\n→ Guardada en {f_comb}")

    # ── Advertencia sobre 18S ─────────────────
    insuf = df18[df18["insuficiente"]][["muestra", "reads", "ambiente"]]
    if not insuf.empty:
        print("\n⚠  Muestras 18S excluidas del resumen por datos insuficientes "
              "(< 10 lecturas tras el procesamiento en mothur):")
        print(insuf.to_string(index=False))


if __name__ == "__main__":
    main()
