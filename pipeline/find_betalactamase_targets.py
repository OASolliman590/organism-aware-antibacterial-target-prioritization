"""Search ChEMBL for beta-lactamase targets and count public activities."""
import requests, time, json
BASE = 'https://www.ebi.ac.uk/chembl/api/data'
queries = ['beta-lactamase', 'KPC-2 beta-lactamase', 'OXA-23 beta-lactamase', 'BlaZ beta-lactamase', 'AmpC beta-lactamase']
found = {}
for q in queries:
    r = requests.get(f'{BASE}/target/search.json', params={'q': q}, timeout=60)
    if r.status_code != 200:
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
    r = requests.get(f'{BASE}/activity', params={'target_chembl_id': tid, 'limit': 1}, timeout=60)
    try:
        d['n_activities'] = r.json()['page_meta']['total_count']
    except Exception:
        d['n_activities'] = 0
    print(tid, d['pref_name'], '|', d['organism'], '|', d['target_type'], '|', d['n_activities'])
    time.sleep(0.25)
json.dump(found, open('/home/ubuntu/work/data/betalactamase_target_hits.json', 'w'), indent=2)
