"""Verify all published interval tables against their expanded-byte identities."""
from pathlib import Path
import hashlib,json
R=Path(__file__).resolve().parents[1]
def main():
    records=json.loads((R/'release/interval_hashes.json').read_text()); n=0
    for p,expected in records.items():
        path=(R/p).resolve()
        if not path.is_relative_to(R):raise ValueError('invalid interval path')
        h=hashlib.sha256()
        with path.open('rb') as f:
            for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
        if h.hexdigest()!=expected:raise ValueError(f'Interval hash mismatch: {p}')
        n+=1
    print(json.dumps({'status':'passed','interval_tables':n}))
if __name__=='__main__':main()
