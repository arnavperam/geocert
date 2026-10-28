from pathlib import Path
import sys,json
from fractions import Fraction as F
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'core'))
from data_api import encode
def save7(name,obj):
 p=ROOT/'run'/name;p.parent.mkdir(exist_ok=True,parents=True);p.write_text(json.dumps(obj,indent=2,default=encode))
