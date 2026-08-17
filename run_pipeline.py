"""Run the complete organism-aware target-prioritization workflow."""
from pathlib import Path
import os
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
os.environ.setdefault('PROJECT_ROOT', str(ROOT))
os.environ.setdefault('INPUT_DIR', str(ROOT / 'inputs'))

steps = [
    ROOT / 'pipeline' / 'prepare_compounds.py',
    ROOT / 'pipeline' / 'target_scoring.py',
    ROOT / 'pipeline' / 'rank_and_figures.py',
    ROOT / 'pipeline' / 'summarize_results.py',
    ROOT / 'pipeline' / 'prior_comparison.py',
]
for step in steps:
    print(f'\n=== {step.name} ===')
    subprocess.run([sys.executable, str(step)], check=True, cwd=ROOT)
print('\nPipeline completed. See results/.')
