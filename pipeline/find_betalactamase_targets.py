"""Refresh a ChEMBL beta-lactamase target diagnostic into a dated snapshot."""
import requests, time, json

try:
    from pipeline.config import load_config
    from pipeline.provenance import utc_now
    from pipeline.snapshots import require_refresh_output, refresh_snapshot_root, write_text_exclusive
except ModuleNotFoundError:
    from config import load_config
    from provenance import utc_now
    from snapshots import require_refresh_output, refresh_snapshot_root, write_text_exclusive


CONFIG = load_config()
CHEMBL_CONFIG = CONFIG.value("refresh.chembl")
OUT = CONFIG.path_for("additional_target_hits")
refresh_snapshot_root(CONFIG)
if str(CHEMBL_CONFIG["source_release"]).lower() == "unrecorded":
    raise RuntimeError(
        "refresh.chembl.source_release must name the ChEMBL release before refresh"
    )
require_refresh_output(CONFIG, OUT)
BASE = 'https://www.ebi.ac.uk/chembl/api/data'
queries = ['beta-lactamase', 'KPC-2 beta-lactamase', 'OXA-23 beta-lactamase', 'BlaZ beta-lactamase', 'AmpC beta-lactamase']
found = {}
gaps = []
for q in queries:
    r = requests.get(f'{BASE}/target/search.json', params={'q': q}, timeout=int(CHEMBL_CONFIG["http_timeout_seconds"]))
    if r.status_code != 200:
        gaps.append({'query': q, 'status': 'request_failed', 'http_status': r.status_code})
        continue
    for t in r.json().get('targets', []):
        tid = t['target_chembl_id']
        if tid not in found:
            found[tid] = {
                'target_chembl_id': tid,
                'pref_name': t.get('pref_name'),
                'organism': t.get('organism'),
                'target_type': t.get('target_type'),
            }
    time.sleep(0.5)
for tid, d in found.items():
    r = requests.get(f'{BASE}/activity', params={'target_chembl_id': tid, 'limit': 1}, timeout=int(CHEMBL_CONFIG["http_timeout_seconds"]))
    try:
        d['n_activities'] = r.json()['page_meta']['total_count']
        d['activity_count_status'] = 'retrieved'
    except Exception:
        d['n_activities'] = None
        d['activity_count_status'] = 'retrieval_failed'
    print(tid, d['pref_name'], '|', d['organism'], '|', d['target_type'], '|', d['n_activities'])
    time.sleep(0.25)
write_text_exclusive(
    CONFIG,
    OUT,
    json.dumps(
        {
            'source': 'ChEMBL',
            'source_release': CHEMBL_CONFIG['source_release'],
            'queried_at_utc': utc_now(),
            'gaps': gaps,
            'targets': found,
        },
        indent=2,
    ) + "\n",
)
