from pathlib import Path
import csv,gzip,json
root=Path('results/proseg_locked')
out={}
for folder in sorted(root.iterdir()):
    if not folder.is_dir():continue
    item={}
    for name in ['transcripts.csv.gz','cells.csv.gz','genes.csv.gz']:
        p=folder/name
        if p.exists():
            with gzip.open(p,'rt') as f:
                rd=csv.DictReader(f);rows=list(rd)
            item[name]={'columns':rd.fieldnames,'rows':len(rows),'preview':rows[:3]}
    if item:out[folder.name]=item
(root/'output_schema.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
