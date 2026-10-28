#!/usr/bin/env python3
"""Retrieve selected ZIP members over validated HTTP ranges; do not fetch images.

This utility prepares external inputs. It is not part of the reported evaluation.
An endpoint without byte-range support fails explicitly instead of silently
fetching a multi-gigabyte archive in response to each seek.
"""
from __future__ import annotations
import argparse,hashlib,io,json,shutil,urllib.request,zipfile
from pathlib import Path

class RemoteZIP(io.RawIOBase):
 def __init__(self,url:str):
  self.url=url;self.pos=0;self.start=0;self.cache=b''
  with urllib.request.urlopen(urllib.request.Request(url,headers={'Range':'bytes=0-0'}),timeout=120) as r:
   if r.status!=206:raise OSError('Server did not honor the initial HTTP Range request')
   cr=r.headers.get('Content-Range','');self.size=int(cr.rsplit('/',1)[1]);self.etag=r.headers.get('ETag');r.read()
 def readable(self):return True
 def seekable(self):return True
 def tell(self):return self.pos
 def seek(self,offset,whence=0):
  v=offset if whence==0 else self.pos+offset if whence==1 else self.size+offset
  if v<0:raise ValueError('negative seek')
  self.pos=v;return v
 def read(self,n=-1):
  if n is None or n<0:n=self.size-self.pos
  n=min(n,max(0,self.size-self.pos))
  if not n:return b''
  if self.start<=self.pos and self.pos+n<=self.start+len(self.cache):
   a=self.pos-self.start;self.pos+=n;return self.cache[a:a+n]
  end=min(self.size,self.pos+max(n,4*1024*1024))-1
  headers={'Range':f'bytes={self.pos}-{end}'}
  if self.etag:headers['If-Match']=self.etag
  with urllib.request.urlopen(urllib.request.Request(self.url,headers=headers),timeout=120) as r:
   expected=f'bytes {self.pos}-{end}/{self.size}'
   if r.status!=206 or r.headers.get('Content-Range')!=expected:raise OSError('Mismatched byte-range response')
   raw=r.read(end-self.pos+2)
  if len(raw)!=end-self.pos+1:raise OSError('Truncated or overlong HTTP range')
  self.start=self.pos;self.cache=raw;self.pos+=n;return raw[:n]

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--url',required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=True)
 stream=RemoteZIP(a.url);records=[]
 wanted={'cells.parquet','transcripts.parquet','cells.zarr.zip','cell_boundaries.parquet','nucleus_boundaries.parquet','gene_panel.json','experiment.xenium'}
 with zipfile.ZipFile(stream) as z:
  members={}
  for info in z.infolist():
   name=Path(info.filename).name
   if name in wanted:
    if name in members:raise ValueError(f'ambiguous member {name}')
    members[name]=info
  for required in ['cells.parquet','transcripts.parquet']:
   if required not in members:raise ValueError(f'missing required member {required}')
  for name,info in sorted(members.items()):
   if info.file_size>5*1024**3:raise ValueError('Member exceeds the declared extraction size guard')
   path=a.out/name;tmp=path.with_suffix(path.suffix+'.partial')
   with z.open(info) as inp,tmp.open('wb') as dst:shutil.copyfileobj(inp,dst,8*1024*1024)
   if tmp.stat().st_size!=info.file_size:raise OSError('decompressed size mismatch')
   digest=hashlib.file_digest(tmp.open('rb'),'sha256').hexdigest();tmp.replace(path)
   records.append({'member':info.filename,'path':name,'bytes':info.file_size,'sha256':digest});print(name,info.file_size,flush=True)
 record={'url':a.url,'archive_bytes':stream.size,'etag':stream.etag,'archive_hash':'not computed: selective member retrieval','members':records}
 (a.out/'download.json').write_text(json.dumps(record,indent=2)+'\n')
if __name__=='__main__':main()
