"""Run pinned Proseg on existing development crops; never call agreement truth.

Install separately: cargo install --locked proseg --version 3.2.0
Requires pandas + pyarrow to read the supplied crop tables.
"""
from __future__ import annotations
import argparse
from collections import Counter
import csv, gzip, hashlib, json, math, shutil, subprocess, time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
VERSION = 'proseg 3.2.0'

def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def load_csv(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', newline='') as f:
        return list(csv.DictReader(f))

def validate_and_compare(source, transcripts, cells):
    """Strict transcript-ID/gene/cell-map join; absent output rows are exclusions."""
    src = {}
    for row in source:
        tid = str(row['transcript_id'])
        if tid in src: raise ValueError('duplicate source transcript ID')
        src[tid] = row
    cmap = {}
    original = set()
    for row in cells:
        key, value = str(row['cell']), str(row['original_cell_id'])
        if key in cmap or value in original: raise ValueError('nonunique cell identity map')
        if not value or value == 'UNASSIGNED': raise ValueError('missing original cell identity')
        cmap[key] = value; original.add(value)
    seen = set(); same = changed = source_unassigned = new_assigned = output_unassigned = 0
    source_counts = Counter(); output_counts = Counter()
    for row in transcripts:
        tid = str(row['transcript_id'])
        if tid in seen: raise ValueError('duplicate output transcript ID')
        if tid not in src: raise ValueError('unknown output transcript ID')
        seen.add(tid); ref = src[tid]
        if row['gene'] != ref['feature_name']: raise ValueError('gene mismatch')
        for axis in ('x', 'y'):
            observed = float(row['observed_' + axis]); expected = float(ref[axis + '_location'])
            if not math.isfinite(observed) or abs(observed-expected) > 0.002:
                raise ValueError('observed coordinate mismatch')
        a = row['assignment']
        if a:
            if a not in cmap: raise ValueError('unknown assigned cell ID')
            owner = cmap[a]
        else:
            owner = 'UNASSIGNED'; output_unassigned += 1
        previous = ref['cell_id']
        same += owner == previous; changed += owner != previous
        source_unassigned += previous == 'UNASSIGNED'
        new_assigned += previous == 'UNASSIGNED' and owner != 'UNASSIGNED'
        if previous != 'UNASSIGNED': source_counts[previous, row['gene']] += 1
        if owner != 'UNASSIGNED': output_counts[owner, row['gene']] += 1
    if not seen: raise ValueError('empty output')
    active = source_counts.keys() | output_counts.keys()
    delta = [abs(output_counts[k]-source_counts[k]) for k in active]
    return {
        'scope': 'descriptive comparison to source assignments on retained output molecules; no truth labels',
        'source_transcripts': len(src), 'output_transcripts': len(seen),
        'excluded_source_ids': sorted(set(src)-seen), 'excluded_source_count': len(src)-len(seen),
        'source_cells_in_retained_rows': len({k[0] for k in source_counts}),
        'proseg_cells': len(cmap), 'same_owner': same, 'changed_owner': changed,
        'source_unassigned_in_retained_rows': source_unassigned,
        'source_unassigned_now_assigned': new_assigned,
        'output_unassigned': output_unassigned,
        'active_cell_gene_pairs': len(active),
        'changed_cell_gene_pairs': sum(d > 0 for d in delta),
        'count_L1_change': sum(delta),
        'accuracy_AUROC_interval_coverage': None,
        'independent_truth_available': False,
    }

def run(binary, crop, out, timeout=900):
    import pandas as pd
    out = Path(out).resolve()
    if out.exists(): raise ValueError('output path exists; use a fresh directory')
    exe = shutil.which(binary)
    if not exe: raise ValueError('Proseg executable unavailable; no experiment run')
    version = subprocess.check_output([exe, '--version'], text=True, timeout=30).strip()
    if version != VERSION: raise ValueError(f'expected {VERSION}, found {version}')
    data = Path(crop).resolve() / 'transcripts.parquet'
    frame = pd.read_parquet(data)
    required = {'transcript_id','cell_id','feature_name','x_location','y_location','z_location','overlaps_nucleus','qv'}
    if not required <= set(frame): raise ValueError('missing Xenium fields')
    if frame.transcript_id.duplicated().any(): raise ValueError('duplicate source transcript ID')
    out.mkdir(parents=True)
    input_path = out / 'input.csv.gz'
    # This derived table is not a new donor or a complete unfiltered Xenium export.
    frame.to_csv(input_path, index=False, compression={'method':'gzip','mtime':0})
    command = [str(Path(exe).resolve()), '--xenium', str(input_path), '--no-factorization',
               '--nthreads', '2', '--output-path', str(out),
               '--output-transcript-metadata','transcript-metadata.csv.gz',
               '--output-cell-metadata','cell-metadata.csv.gz',
               '--output-gene-metadata','gene-metadata.csv.gz',
               '--output-counts','counts.mtx.gz']
    record = {'status':'running', 'version':version, 'binary_sha256':digest(exe),
              'input_parquet_sha256':digest(data), 'input_csv_sha256':digest(input_path),
              'command':command, 'threads':2, 'factorization':False,
              'sampling':'upstream defaults: burnin 200, samples 200, hillclimb 50; stochastic, no seed option',
              'model_change_reason':'default factorization failed on roi_01 small crop with EmptyInput; adaptation disclosed before comparison',
              'data_role':'existing development donor; crop edges and prefiltering limit interpretation'}
    path = out / 'run.json'; path.write_text(json.dumps(record,indent=2)+'\n')
    start=time.monotonic()
    try:
        with (out/'run.log').open('w') as log:
            proc=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,timeout=timeout)
        record.update(returncode=proc.returncode, elapsed_seconds=time.monotonic()-start)
        if proc.returncode: raise RuntimeError('Proseg exited unsuccessfully; inspect run.log')
        source = load_csv(input_path)
        record['comparison'] = validate_and_compare(source,load_csv(out/'transcript-metadata.csv.gz'),load_csv(out/'cell-metadata.csv.gz'))
        record['output_sha256']={p.name:digest(p) for p in out.iterdir() if p.is_file() and p.name not in {'run.json','run.log'}}
        record['status']='completed_and_validated'
    except Exception as exc:
        record.update(status='failed',error=str(exc),elapsed_seconds=time.monotonic()-start)
        raise
    finally:
        path.write_text(json.dumps(record,indent=2)+'\n')
    return record

if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary',default='proseg')
    parser.add_argument('--crop',required=True,type=Path)
    parser.add_argument('--output',required=True,type=Path)
    parser.add_argument('--timeout',type=int,default=900)
    a=parser.parse_args()
    print(json.dumps(run(a.binary,a.crop,a.output,a.timeout),indent=2))
