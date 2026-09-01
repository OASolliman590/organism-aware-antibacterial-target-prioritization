"""TEMPORARY diagnostic: locate where the v2 golden differs between hosts.

Prints per-column digests and reference counts so two machines can be compared
directly. Delete once the v2 golden reproduces across hosts.
"""

from __future__ import annotations

import hashlib
import platform
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from pipeline import open_target_discovery_v2 as v2  # noqa: E402
from tests.test_v2_golden import GOLDEN_DECIMALS, _public_v2_outputs  # noqa: E402


def digest(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()[:16]


def main() -> None:
    print("platform:", platform.platform(), "| python:", sys.version.split()[0])
    import numpy, rdkit  # noqa: E402

    print("numpy:", numpy.__version__, "| rdkit:", rdkit.__version__)
    print("numpy BLAS:", (numpy.__config__.show(mode="dicts") or {}).get("Build Dependencies", {}).get("blas", {}).get("name", "?"))

    refs = v2.load_refs()
    print("\n--- reference counts per class ---")
    for key in sorted(refs):
        print(f"{key:24s} {len(refs[key])}")

    frame = _public_v2_outputs()["v2_benchmark_open_target_scores.csv"]
    rounded = frame.round(GOLDEN_DECIMALS)
    print("\n--- per-column digest (rounded) ---")
    for column in rounded.columns:
        print(f"{column:38s} {rounded[column].to_csv(index=False)[:0] or ''}"
              f"{digest(rounded[column].to_csv(index=False))}")

    print("\n--- first 3 rows ---")
    with pd.option_context("display.max_columns", None, "display.width", 300):
        print(rounded.head(3).to_string())


if __name__ == "__main__":
    main()
