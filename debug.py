from pathlib import Path

from polysaccharide.polymers.management import PolymerManager
from polysaccharide.polymers.building import build_linear_polymer, mbmol_from_mono_smarts

p = Path('Collections')

mgr = PolymerManager(p / 'water_soluble_large')
pdir = mgr.polymers['peg_modified']
mono_smarts = pdir.monomer_info.monomers

for name, smarts in mono_smarts.items():
    print(name)
    mbmol = mbmol_from_mono_smarts(smarts)

chain = build_linear_polymer(monomer_smarts=mono_smarts, DOP=10)
print(chain)