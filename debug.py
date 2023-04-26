from polymer_utils.charging.averaging import write_lib_chgs_from_mono_data
from polymer_utils.representation import PolymerManager
from pathlib import Path

main_ff_xml = Path('resources')/'force_fields'/'openff_constrained-2.0.0.offxml'
p = Path('Collections') / 'simple_polymers_updated'
mgr = PolymerManager(p)
pdir = mgr.polymers['naturalrubber_solv_water']

ff, lc = write_lib_chgs_from_mono_data(pdir.monomer_data_charged, main_ff_xml, Path('test.offxml'))