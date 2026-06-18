# src/gaussflow/cli.py

from __future__ import annotations

import argparse
from pathlib import Path

from gaussflow.workflow import submit_workflow


def main():
    parser = argparse.ArgumentParser(
        prog="gaussflow",
        description="Minimal Gaussian workflow submitter.",
    )

    subparsers = parser.add_subparsers(dest="command", required=True)

    submit_parser = subparsers.add_parser(
        "submit",
        help="Generate Gaussian input files and submit jobs.",
    )
    submit_parser.add_argument(
        "--config",
        required=True,
        type=Path,
        help="Path to config JSON file.",
    )

    args = parser.parse_args()

    if args.command == "submit":
        submit_workflow(config_path=args.config)


if __name__ == "__main__":
    main()