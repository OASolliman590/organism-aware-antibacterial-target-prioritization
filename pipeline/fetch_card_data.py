"""Download and extract the current public CARD archives when missing.

The raw CARD release is deliberately not committed to Git because it is large.
Derived v2 resistance tables may be retained in a local analysis bundle; a clean
public checkout can regenerate them with this script when internet access exists.
"""
from pathlib import Path
import bz2, io, tarfile, urllib.request

try:
    from pipeline.config import load_config
    from pipeline.snapshots import refresh_snapshot_root, require_refresh_output, write_bytes_exclusive
except ModuleNotFoundError:  # direct ``python pipeline/<script>.py`` execution
    from config import load_config
    from snapshots import refresh_snapshot_root, require_refresh_output, write_bytes_exclusive


CONFIG = load_config()
BASE = CONFIG.path_for("card_raw")
DATA=BASE/'data'; ONTO=BASE/'ontology'
CARD_CONFIG = CONFIG.value("refresh.card")
URLS={'card-data.tar.bz2':'https://card.mcmaster.ca/latest/data','card-ontology.tar.bz2':'https://card.mcmaster.ca/latest/ontology'}

def fetch(name,url):
    out=BASE/name
    require_refresh_output(CONFIG, out)
    print('downloading',url,flush=True)
    req=urllib.request.Request(url,headers={'User-Agent':'v3-antibacterial-target-prioritization/1.0'})
    with urllib.request.urlopen(req,timeout=int(CARD_CONFIG["request_timeout_seconds"])) as response:
        payload = response.read()
    if len(payload) < 1000:
        raise RuntimeError(f"CARD download was unexpectedly small: {url}")
    write_bytes_exclusive(CONFIG, out, payload)
    return out

def extract(path):
    dest=ONTO if 'ontology' in path.name else DATA
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
                    write_bytes_exclusive(CONFIG, target, src.read())

def main():
    refresh_snapshot_root(CONFIG)
    if str(CARD_CONFIG["source_version"]).lower() == "unrecorded":
        raise RuntimeError(
            "refresh.card.source_version must name the CARD release before refresh"
        )
    for name,url in URLS.items(): extract(fetch(name,url))
    print('CARD data ready:',DATA,'ontology:',ONTO)

if __name__=='__main__': main()
