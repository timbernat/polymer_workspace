from pathlib import Path
from shutil import copyfile

from polymer_utils.representation import PolymerManager
from polymer_utils.molutils.building import build_linear_polymer_limited

COLL_PATH = Path('Collections')
COMPAT_PDB_PATH = Path('compatible_pdbs_updated')

# Set parameters here
chain_lim = 180
source_coll_name = 'water_soluble_polymers'
output_coll_name = 'water_soluble_small'
flip_term_group_labels = ['paam_modified']

if __name__ == '__main__':
    # load source polymers
    poly_source_path = COMPAT_PDB_PATH / source_coll_name
    collection_path  = COLL_PATH / source_coll_name

    structure_path   = poly_source_path / f'{poly_source_path.name}_structures'
    monomer_path     = poly_source_path / f'{poly_source_path.name}_monomers'

    mgr = PolymerManager(collection_path)
    if not mgr.mol_dirs: 
        mgr.populate_collection(struct_dir=structure_path, monomer_dir=monomer_path) # ensure originals have been loaded
    
    # create output dirs
    reduced_dir = COMPAT_PDB_PATH / output_coll_name
    reduced_structures = reduced_dir / f'{reduced_dir.name}_structures'
    reduced_monomers   = reduced_dir / f'{reduced_dir.name}_monomers'

    reduced_dir.mkdir(       exist_ok=True)
    reduced_monomers.mkdir(  exist_ok=True)
    reduced_structures.mkdir(exist_ok=True)

    # generate new structures
    reverse = False # needed
    for pdir in mgr.mol_dirs_list:
        print(pdir.mol_name)
        monomers = pdir.monomer_data['monomers']
        if pdir.mol_name in flip_term_group_labels:
            monomers.pop('paam_SPECIAL_TERM')   # special case needed for extra messy term group in PAAM...
            reverse = True                      # ...TODO : find generalized way to isolate head and tail groups from larger collection of monomers

        chain = build_linear_polymer_limited(monomers, max_chain_len=chain_lim, reverse_term_labels=reverse)
        chain.save(str(reduced_structures/f'{pdir.mol_name}.pdb'), overwrite=True)
        copyfile(pdir.monomer_file, reduced_monomers/f'{pdir.mol_name}.json')