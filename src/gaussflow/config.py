# src/gaussflow/config.py

from __future__ import annotations

import json
from pathlib import Path
from typing import Any


def load_config(config_path):
    config_path = Path(config_path)

    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")

    with config_path.open("r", encoding="utf-8") as f:
        config = json.load(f)

    if not isinstance(config, dict):
        raise ValueError("Config must be a JSON object.")

    return config

def dump_config(config, config_path):
    config_path = Path(config_path)

    with config_path.open("w", encoding="utf-8") as f:
        json.dump(config, f)
