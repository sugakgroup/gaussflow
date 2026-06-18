
import argparse
import csv
from pathlib import Path
import json

from gaussflow.molecule import load_molecules


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=str, required=True)

    config_file = parser.parse_args().config
    config = json.loads(Path(config_file).read_text(encoding="utf-8"))

    molecules = load_molecules(config)

    with open(Path(f"{config['output']['root']}/summary.csv"), "w", encoding="utf-8") as f_out:
        writer = csv.writer(f_out)
        writer.writerow([col for col in config["after_jobs"]["summary_columns"]])
        for mol in molecules:
            mol_id = mol["id"]
            row = []
            with open(Path(f"{config['output']['root']}/mol_{mol_id}/properties.json"), "r", encoding="utf-8") as f_in:
                properties = json.load(f_in)
                for col in config["after_jobs"]["summary_columns"]:
                    item = properties.copy()
                    for key in col.split(":")[::-1]:
                        item = item[key]
                    row.append(str(item))
            writer.writerow(row)

if __name__ == "__main__":
    main()