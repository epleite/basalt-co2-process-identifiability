from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
STEPS = [
    ROOT / 'code' / 'reproduce_v1_3_fem.py',
    ROOT / 'code' / 'reproduce_v1_4.py',
]

for step in STEPS:
    print(f'\n=== {step.name} ===', flush=True)
    subprocess.run([sys.executable, str(step)], cwd=ROOT, check=True)

print('\nExecutable GitHub computational core completed successfully.')
