"""Count compatibility with a published segmenter, not biological accuracy."""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter,defaultdict
import json,gzip,csv,sys
from shapely.geometry import shape,Point
from shapely.strtree import STRtree
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from global_model.inference import BUDGETS

def readcsv(p):
    with gzip.open(p,'rt') as f:return list(csv.DictReader(f))
def readjson(p):
    with gzip.open(p,'rt') as f:return json.load(f)

def run():
    source=json.loads((ROOT/'data/observed_crop_counts.json').read_text())['rois']
    spectra=json.loads((ROOT/'results/cell_spectra.json').read_text());index={(r['roi'],r['cell'],r['gene']):r for r in spectra}
    family=json.loads((ROOT/'data/global_patches.json').read_text());membership={}
    for p in family['patches']:
        for z in p['points']:membership[p['roi'],z['transcript_id']]=p
    complete={r['roi']:{c['cell'] for c in r['cells'] if c['complete_observed_cell']} for r in family['census']}
    results=[];rows=[]
    for setting in ['proseg_locked','proseg_unfactorized']:
        base=ROOT/'results'/setting
        if not base.exists():continue
        for folder in sorted(base.iterdir()):
            if not folder.is_dir() or not (folder/'transcripts.csv.gz').exists():continue
            roi=next((r for r in source if folder.name.startswith(r+'_')),None)
            if roi is None:continue
            meta=json.loads((ROOT/'data/source_crops'/roi/'metadata.json').read_text())
            name_to_label={v:int(k) for k,v in meta['label_to_cell_id'].items()}
            cells=readcsv(folder/'cells.csv.gz');out_to_raw={int(c['cell']):name_to_label[c['original_cell_id']] for c in cells if c['original_cell_id'] in name_to_label}
            missing=sorted(complete[roi]-set(out_to_raw.values()))
            transcripts=readcsv(folder/'transcripts.csv.gz');byid={z['transcript_id']:z for z in transcripts}
            counts_reported=defaultdict(Counter);owner_reported={};unknown=0
            for z in transcripts:
                owner=0 if not z['assignment'] or z['background'].lower()=='true' else out_to_raw.get(int(z['assignment']),-1)
                if owner<0:unknown+=1;continue
                owner_reported[z['transcript_id']]=owner
                if owner:counts_reported[owner][z['gene']]+=1
            polys=readjson(folder/'polygons.geojson.gz');features=polys['features']
            geoms=[];labels=[]
            for z in features:
                pr=z.get('properties',{});raw=pr.get('cell',pr.get('cell_id',None))
                if raw is None:raise ValueError('Unknown Proseg polygon cell key: '+str(pr))
                k=out_to_raw.get(int(raw))
                if k is None:continue
                geom=shape(z['geometry'])
                if not geom.is_valid:raise ValueError('Invalid Proseg polygon')
                geoms.append(geom);labels.append(k)
            tree=STRtree(geoms);counts_poly=defaultdict(Counter);owner_poly={};overlap=0
            for z in source[roi]['records']:
                pt=Point(float(z['x_um']),float(z['y_um']));cand=tree.query(pt)
                found=[labels[int(i)] for i in cand if geoms[int(i)].covers(pt)]
                if len(set(found))>1:overlap+=1;continue
                owner=found[0] if found else 0;owner_poly[z['id']]=owner
                if owner:counts_poly[owner][z['gene']]+=1
            for mode,cc,oo in [('reported_assignment',counts_reported,owner_reported),('polygon_at_original_coordinates',counts_poly,owner_poly)]:
                rr=[];changes=0;outside=0;unsupported=0;unobserved=0
                for z in source[roi]['records']:
                    if z['id'] not in oo:unobserved+=1;continue
                    dest=oo[z['id']]
                    if dest!=z['cell']:
                        changes+=1;p=membership.get((roi,z['id']))
                        if p is None:outside+=1
                        elif dest not in p['cell_labels']:unsupported+=1
                for cell in sorted(complete[roi]&set(out_to_raw.values())):
                    genes=set(source[roi]['counts'].get(str(cell),{}))|set(cc[cell])|{g for r,k,g in index if r==roi and k==cell}
                    for gene in sorted(genes):
                        target=cc[cell].get(gene,0);s=index.get((roi,cell,gene))
                        spec={int(k):F(v) for k,v in s['minima'].items()} if s else {0:F(0)}
                        lower=min(spec);upper=max(spec);exactcost=spec.get(target)
                        hulldelta=min((max(a,b) for k,a in spec.items() for l,b in spec.items() if k<=target<=l),default=None)
                        r={'method_run':setting+'/'+folder.name,'readout':mode,'roi':roi,'cell':cell,'gene':gene,'reference_count':source[roi]['counts'].get(str(cell),{}).get(gene,0),
                           'method_count':target,'exact_count_min_area':None if exactcost is None else str(exactcost),
                           'interval_containment_min_area':None if hulldelta is None else str(hulldelta),'unbudgeted_lower':lower,'unbudgeted_upper':upper}
                        rr.append(r);rows.append(r)
                summary={'method_run':setting+'/'+folder.name,'readout':mode,'questions':len(rr),'complete_cells_compared':len(complete[roi]&set(out_to_raw.values())),
                         'complete_cells_missing_identity':missing,'input_transcripts':len(source[roi]['records']),'output_transcripts':len(transcripts),
                         'input_records_absent_from_assignment_table':len(set(z['id'] for z in source[roi]['records'])-set(byid)),
                         'unknown_output_labels':unknown,'polygon_overlap_records':overlap,
                         'changed_observed_owners':changes,'changes_outside_moving_patches':outside,'changes_to_nonpatch_label':unsupported,'unobserved_in_readout':unobserved,
                         'exact_counts_unattainable':sum(r['exact_count_min_area'] is None for r in rr),
                         'count_outside_all_budget_intervals':sum(r['interval_containment_min_area'] is None for r in rr),'budgets':[]}
                for d in BUDGETS:
                    summary['budgets'].append({'budget':str(d),'interval_contains':sum(r['interval_containment_min_area'] is not None and F(r['interval_containment_min_area'])<=d for r in rr),
                                               'exact_count_attainable':sum(r['exact_count_min_area'] is not None and F(r['exact_count_min_area'])<=d for r in rr)})
                results.append(summary)
    output={'scope':'Retained observed counts compared to stochastic Proseg point estimates and consensus polygons. Not ground truth, coverage probability, or parameter calibration.','runs':results}
    (ROOT/'results/proseg_comparison_summary.json').write_text(json.dumps(output,indent=2)+'\n')
    with (ROOT/'results/proseg_count_comparison.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    print(json.dumps(output,indent=2))
if __name__=='__main__':run()
