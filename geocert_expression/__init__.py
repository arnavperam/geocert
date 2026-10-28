"""Model-conditioned expression certificates, candidate release 0.2.0a1."""
from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[1]
ENGINE = Path(__file__).resolve().parent / "_engine"
sys.path.insert(0, str(ENGINE))
__version__ = "0.2.0a1"
