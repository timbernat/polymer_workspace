'''For transferring charged monomer information to full-sized WaSPs once ABE10 charging is done on reduced WaSPs'''

# Generic imports
from pathlib import Path

# Custom Imports
from polysaccharide.representation import PolymerManager
from polysaccharide.solvation.solvents import WATER_TIP3P

# Static Paths
RESOURCE_PATH = Path('resources')
COLL_PATH = Path('Collections')
COMPAT_PDB_PATH = Path('compatible_pdbs_updated')

# ------------------------------------------------------------------------------

# Set parameters here
small_coll = COLL_PATH / 'water_soluble_small'
new_coll = COLL_PATH / 'water_soluble_large'
structure_dir = COMPAT_PDB_PATH / 'water_soluble_polymers' / 'water_soluble_polymers_structures' 

if __name__ == '__main__':
    new_coll.mkdir(exist_ok=True)

    # load small collection and clone into new folder
    mgr_small = PolymerManager(small_coll)
    for polymer in mgr_small.polymers_list:
        if polymer.solvent == WATER_TIP3P:
            polymer.clone(
                dest_dir=new_coll,
                clone_name=polymer.base_mol_name,
                clone_solvent=True,
                clone_structures= False,
                clone_monomers=True,
                clone_ff= False,
                clone_charges= False,
                clone_sims= False
            )

    # copy large structures over to clones, initialize as PolymerManager
    mgr_large = PolymerManager(new_coll)
    for polymer in mgr_large.polymers_list: # copy large structuers over
        polymer.populate_pdb(structure_dir)