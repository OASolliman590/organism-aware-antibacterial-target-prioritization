"""Download and extract the current public CARD archives when missing.

The raw CARD release is deliberately not committed to Git because it is large.
Derived v2 resistance tables may be retained in a local analysis bundle; a clean
public checkout can regenerate them with this script when internet access exists.
"""
from pathlib import Path
import bz2, io, os, tarfile, urllib.request

ROOT=Path(os.environ.get('PROJECT_ROOT',Path(__file__).resolve().parents[1])); BASE=ROOT/'data'/'card'; DATA=BASE/'data'; ONTO=BASE/'ontology'
URLS={'card-data.tar.bz2':'https://card.mcmaster.ca/latest/data','card-ontology.tar.bz2':'https://card.mcmaster.ca/latest/ontology'}

def fetch(name,url):
    out=BASE/name
    if not out.exists() or out.stat().st_size<1000:
        print('downloading',url,flush=True)
        req=urllib.request.Request(url,headers={'User-Agent':'v2-antibacterial-target-discovery/1.0'})
        with urllib.request.urlopen(req,timeout=60) as r, open(out,'wb') as f:
            while True:
                chunk=r.read(1024*1024)
                if not chunk: break
                f.write(chunk)
    return out

def extract(path):
    dest=ONTO if 'ontology' in path.name else DATA
    dest.mkdir(parents=True,exist_ok=True)
    with bz2.open(path,'rb') as stream:
        raw=stream.read()
    with tarfile.open(fileobj=io.BytesIO(raw),mode='r:') as tar:
        for member in tar.getmembers():
            name=Path(member.name).name
            if not name or name in {'.','..'}: continue
            # The current archives contain flat files in their respective folders.
            target=dest/name
            if member.isfile():
                src=tar.extractfile(member)
                if src:
                    target.write_bytes(src.read())

def main():
    BASE.mkdir(parents=True,exist_ok=True)
    for name,url in URLS.items(): extract(fetch(name,url))
    print('CARD data ready:',DATA,'ontology:',ONTO)

if __name__=='__main__': main()
