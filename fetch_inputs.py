#!/usr/bin/env python3
"""Fetch checksum-pinned public inputs without overwriting scientific source.

The input archive is immutable by commit and SHA-256. It contains the same lung
and breast format-test crops as this study, plus the imported Proseg executions.
Only data and baseline result files are installed; source programs are excluded.
"""
from pathlib import Path, PurePosixPath
import argparse, hashlib, io, json, os, shutil, stat, tempfile, urllib.request, zipfile

ROOT=Path(__file__).resolve().parent
REVISION='2f42eef0afab163502470c77b339f5b8401d1add'
BASE=f'https://raw.githubusercontent.com/arnavperam/geocert/{REVISION}/atlas2027/artifact_parts/'
ARCHIVE_SHA='749e55845d9f5d47f2b0f4ff17bda3a0cd33b7c40775223dd32e42d1f813f383'
ARCHIVE_BYTES=24431242
PARTS=[('part000',4194304,'2ffe0e09202c96c18d2e699989530dc019736d8e8bc64ce1ec6ef1b13988d2af'),
 ('part001',4194304,'e6e8aefe6219100dd3059f214d2286961d8a010db6eae766f2ccceeb43d5df3e'),
 ('part002',4194304,'70b070c2841dc8b86da28b94607c70d7a7bfd2e61e6e746745f668c87db137c7'),
 ('part003',4194304,'d18badbdfdbbd824c796085a8295d34413ec8dd33100bcbc9dac719767422bdf'),
 ('part004',4194304,'5cda66249554e781da5355badab9fd018705420832f15295b3b80d185f37b2dc'),
 ('part005',3459722,'6a5ddb3b83467c42af4af649d345f10d4ef091cf5347fc25628e0175a172e6b8')]

def digest(b):return hashlib.sha256(b).hexdigest()

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--archive',type=Path,help='Use an already downloaded exact archive');args=ap.parse_args()
 if args.archive:raw=args.archive.read_bytes()
 else:
  pieces=[]
  for name,size,sha in PARTS:
   req=urllib.request.Request(BASE+name,headers={'User-Agent':'GeoCert-reproduction/1.1.0'})
   with urllib.request.urlopen(req,timeout=120) as f:b=f.read(size+1)
   if len(b)!=size or digest(b)!=sha:raise ValueError('input part identity mismatch: '+name)
   pieces.append(b)
  raw=b''.join(pieces)
 if len(raw)!=ARCHIVE_BYTES or digest(raw)!=ARCHIVE_SHA:raise ValueError('input archive identity mismatch')
 installed=[]
 with zipfile.ZipFile(io.BytesIO(raw)) as z:
  if z.testzip() is not None:raise ValueError('input archive CRC failed')
  for entry in z.infolist():
   p=PurePosixPath(entry.filename)
   if p.is_absolute() or '..' in p.parts:raise ValueError('unsafe archive path')
   if entry.is_dir():continue
   if len(p.parts)<2:continue
   q=PurePosixPath(*p.parts[1:]);mode=entry.external_attr>>16
   if stat.S_ISLNK(mode):raise ValueError('input archive symlink is unsupported')
   selected=q.parts[0]=='data' or q.parts[:2] in [('external','data'),('external','results')]
   if q.parts[0]=='results' and len(q.parts)>1:
    selected=q.parts[1] not in ['revision_reach','proseg']
   if q.parts[:2]==('results','proseg'):
    selected=True;q=PurePosixPath('data','proseg',*q.parts[2:])
   if not selected:continue
   target=ROOT/Path(q);target.parent.mkdir(parents=True,exist_ok=True);content=z.read(entry)
   if target.exists():
    if target.read_bytes()!=content:
     # Existing scientifically valid local records are preserved, never overwritten.
     continue
   else:
    fd,tmp=tempfile.mkstemp(prefix='geocert-input-',dir=target.parent)
    try:
     with os.fdopen(fd,'wb') as f:f.write(content)
     os.replace(tmp,target)
    finally:
     if os.path.exists(tmp):os.unlink(tmp)
   installed.append(str(q))
 expected=ROOT/'docs/BASELINE_INPUT_HASHES.json'
 if expected.exists():
  for relative,sha in json.loads(expected.read_text()).items():
   p=ROOT/relative
   if not p.is_file() or digest(p.read_bytes())!=sha:raise ValueError('prepared input hash mismatch: '+relative)
 print(json.dumps({'archive_sha256':ARCHIVE_SHA,'input_files_present':len(installed),'source_files_overwritten':0,'status':'passed'}))
if __name__=='__main__':main()
