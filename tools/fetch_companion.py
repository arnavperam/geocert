"""Download and verify the unchanged full publication companion archive.

Usage: python tools/fetch_companion.py [output.zip]
The two stored parts avoid an upload transport limit; concatenation reproduces
the original ZIP exactly. No archive is extracted or executed by this command.
"""
from pathlib import Path
import hashlib,sys,urllib.request

EXPECTED='c786d48ce8925c5e5fe45dd4a96a7782edc5ba61377d2f833e9ddd0c9c53ba60'
REF='reader-revision'
BASE=f'https://raw.githubusercontent.com/arnavperam/geocert/{REF}/artifacts/'
PARTS=['publication_companion.zip.part01','publication_companion.zip.part02']

def main():
    target=Path(sys.argv[1]) if len(sys.argv)>1 else Path('run/publication_companion.zip')
    if target.exists():
        if hashlib.sha256(target.read_bytes()).hexdigest()==EXPECTED:
            print('Already verified:',target);return
        raise SystemExit('Refusing to overwrite a different existing file: '+str(target))
    target.parent.mkdir(parents=True,exist_ok=True);temp=target.with_name(target.name+'.partial')
    if temp.exists():raise SystemExit('Remove or rename the existing partial file before retrying: '+str(temp))
    h=hashlib.sha256()
    with temp.open('xb') as out:
        for name in PARTS:
            with urllib.request.urlopen(BASE+name,timeout=60) as response:
                while chunk:=response.read(1024*1024):h.update(chunk);out.write(chunk)
    if h.hexdigest()!=EXPECTED:raise SystemExit('Companion checksum mismatch; partial file retained for inspection.')
    temp.rename(target);print('Verified full companion:',target,'SHA256',EXPECTED)

if __name__=='__main__':main()
