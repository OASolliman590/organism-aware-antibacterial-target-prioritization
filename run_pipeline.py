"""Run the v2 open-target-discovery workflow.

Private structures are optional and are read only from local ignored paths.
Public benchmark/reference/annotation modules remain reproducible without private
inputs when their public data are available.
"""
from pathlib import Path
import os
import subprocess
import sys

ROOT = Path(os.environ.get('PROJECT_ROOT', Path(__file__).resolve().parent))
os.environ.setdefault('PROJECT_ROOT', str(ROOT))
os.environ.setdefault('INPUT_DIR', str(ROOT / 'inputs'))

def run(rel):
    step=ROOT/rel
    print(f'\n=== {step.name} ===',flush=True)
    subprocess.run([sys.executable,str(step)],check=True,cwd=ROOT,env=os.environ.copy())

steps=[
    'pipeline/fetch_benchmark_structures.py',
    'pipeline/fetch_chembl_reference_subtypes_v21.py',
    'pipeline/build_reference_quality.py',
    'pipeline/fetch_card_data.py',
    'pipeline/fetch_species_targets.py',
    'pipeline/sequence_compatibility.py',
    'pipeline/build_card_resistance_annotations.py',
    'pipeline/parse_card_snps_v2.py',
    'pipeline/fetch_structure_catalog_v2.py',
    'pipeline/open_target_discovery_v2.py',
    'pipeline/benchmark_v2.py',
    'pipeline/calibrate_uncertainty_v2.py',
    'pipeline/build_validation_plan_v2.py',
    'pipeline/v2_figures.py',
    'pipeline/summarize_v2.py',
]
for rel in steps:
    if (ROOT/rel).exists(): run(rel)
print('\nV2.1 open-target-discovery pipeline completed. See results/ for local outputs and data/ for public annotations.')
