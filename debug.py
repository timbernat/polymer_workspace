from polysaccharide.charging.averaging import write_lib_chgs_from_mono_data
from polysaccharide.representation import PolymerManager
from pathlib import Path

main_ff_xml = Path('resources')/'force_fields'/'openff_constrained-2.0.0.offxml'
p = Path('Collections') / 'water_soluble_large'
mgr = PolymerManager(p)
pdir = mgr.polymers['pnipam_modified_solv_water']

print(pdir.offmol)