""" EMOD′ with the fixes named by the importing file (variants/<fix>.py; `all.py` every fix): its factories re-exported """
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / 'models' / 'emod' / 'v4'))
