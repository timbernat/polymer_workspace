from pathlib import Path

from polysaccharide.representation import PolymerManager
from polysaccharide.molutils.polymer.building import build_linear_polymer, mbmol_from_mono_smarts

p = Path('Collections')

mgr = PolymerManager(p / 'water_soluble_large')
pdir = mgr.polymers['peg_modified']
mono_smarts = pdir.monomer_data['monomers']

for name, smarts in mono_smarts.items():
    print(name)
    mbmol = mbmol_from_mono_smarts(smarts)

chain = build_linear_polymer(monomer_smarts=mono_smarts, DOP=10)
print(chain)