import argparse
import sys
from pathlib import Path

METAL_RESNAMES = {
    "AL",
    "BA",
    "CA",
    "CD",
    "CO",
    "CS",
    "CU",
    "FE",
    "HG",
    "K",
    "MG",
    "MN",
    "NA",
    "NI",
    "PB",
    "SR",
    "ZN",
}

METAL_ELEMENTS = {
    "AL",
    "BA",
    "CA",
    "CD",
    "CO",
    "CS",
    "CU",
    "FE",
    "HG",
    "K",
    "MG",
    "MN",
    "NA",
    "NI",
    "PB",
    "SR",
    "ZN",
}


def is_metal_line(line: str) -> bool:
    if not (line.startswith("HETATM") or line.startswith("ATOM")):
        return False
    resname = line[17:20].strip().upper()
    element = line[76:78].strip().upper()
    if resname in METAL_RESNAMES:
        return True
    if element in METAL_ELEMENTS:
        return True
    return False


def preprocess_receptor(input_path: Path, output_path: Path) -> tuple[int, int]:
    removed = 0
    total = 0
    with input_path.open("r", encoding="utf-8", errors="replace") as infile, output_path.open(
        "w", encoding="utf-8"
    ) as outfile:
        for line in infile:
            total += 1
            if is_metal_line(line):
                removed += 1
                continue
            outfile.write(line)
    return removed, total


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Preprocesa un PDB eliminando iones/metales comunes."
    )
    parser.add_argument("input_pdb", type=Path, help="PDB de entrada")
    parser.add_argument("output_pdb", type=Path, help="PDB de salida preprocesado")
    args = parser.parse_args()

    removed, total = preprocess_receptor(args.input_pdb, args.output_pdb)
    if removed:
        print(
            f"Se eliminaron {removed} líneas con iones/metales de {total} líneas.",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
