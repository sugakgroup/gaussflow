# src/gaussflow/molecule.py

from pathlib import Path
import csv

def load_molecules(config):
    input_config = config["input"]

    table_path = Path(input_config["molecule_table"])

    id_col = input_config.get("id_column")

    molecules = []

    with table_path.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)

        for i, row in enumerate(reader):
            molecules.append({key: value for key, value in row.items() if key != id_col})
            if id_col is None and "id" not in molecules[-1]:
                molecules[-1]["id"] = i

    return molecules