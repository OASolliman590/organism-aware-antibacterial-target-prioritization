"""Run the revised open target-discovery workflow.

Private structures are optional and are read only from the local ignored
inputs/data/compounds directories. Public benchmark/reference data are processed
regardless of whether private inputs are present.
"""
from pathlib import Path
import os
import subprocess
import sys

ROOT = Path(os.environ.get('PROJECT_ROOT', Path(__file__).resolve().parent))
os.environ.setdefault('PROJECT_ROOT', str(ROOT))
os.environ.setdefault('INPUT_DIR', str(ROOT / 'inputs'))

steps = [
    ROOT / 'pipeline' / 'fetch_benchmark_structures.py',
    ROOT / 'pipeline' / 'open_target_discovery.py',
    ROOT / 'pipeline' / 'summarize_benchmark.py',
    ROOT / 'pipeline' / 'open_target_figures.py',
]
for step in steps:
    print(f'\n=== {step.name} ===')
    subprocess.run([sys.executable, str(step)], check=True, cwd=ROOT, env=os.environ.copy())
print('\nOpen target-discovery pipeline completed. See results/.')
