from pathlib import Path
import os

ROOT=Path(os.environ.get('PROJECT_ROOT',Path(__file__).resolve().parents[1]))
patterns=['README.md','run_pipeline.py','.gitignore','docs/**/*.md','data/**/*.csv','data/**/*.json','data/**/*.fasta','pipeline/*.py']
paths=[]
for pattern in patterns: paths.extend(ROOT.glob(pattern))
for p in sorted(set(paths)):
    if not p.is_file(): continue
    raw=p.read_bytes()
    if b'\x00' in raw: continue
    clean=raw.replace(b'\r\n',b'\n').replace(b'\r',b'\n')
    lines=clean.split(b'\n')
    clean=b'\n'.join(line.rstrip(b' \t') for line in lines)
    if clean!=raw: p.write_bytes(clean)
print('sanitized',len(set(paths)),'public text candidates')
