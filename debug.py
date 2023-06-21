from pathlib import Path

from polysaccharide.polymer.management import PolymerManager

p = Path('Collections')

mgr = PolymerManager(p / 'water_soluble_polymers_equil')
pdir = mgr.polymers['paam_modified_conf_1_solv_water_equil']
offmol = pdir.offmol