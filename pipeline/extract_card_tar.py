from pathlib import Path
import bz2, tarfile, sys

for archive, outdir in [(Path('data/card/card-ontology.tar.bz2'),Path('data/card/ontology')),(Path('data/card/card-data.tar.bz2'),Path('data/card/data'))]:
    if not archive.exists():
        print('missing',archive); continue
    outdir.mkdir(parents=True,exist_ok=True)
    try:
        with bz2.open(archive,'rb') as fh:
            with tarfile.open(fileobj=fh,mode='r|') as tar:
                tar.extractall(outdir)
        print('extracted',archive,'to',outdir)
    except Exception as exc:
        print('failed',archive,repr(exc))
        sys.exit(1)
