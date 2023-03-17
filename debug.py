from pathlib import Path
from polymer_utils.representation import PolymerDir
from openff.toolkit.topology import Topology


def load_mol_and_topo(pdb_path : Path, json_path : Path, verbose : bool=False):
    '''Load Molecule and Topology from a pdb and a monomer json file, performing residue matching on monomer units
    Assumes that the pdb only contains has ONE simple homopolymer (will only load first molecule if multiple are present'''
    off_topology, _, error = Topology.from_pdb_and_monomer_info(str(pdb_path), json_path, strict=True, verbose=verbose)
    mol = next(off_topology.molecules) # get the first molecule (assumed to be the polymer of interest)

    return mol, off_topology


solvent = 'water'

solv_src = Path('Core/solvents')/solvent
solvent_pdb = solv_src/f'{solvent}.pdb'
solvent_json = solv_src/f'{solvent}.json'



pdir = PolymerDir.from_file(Path('Polymers/polyvinylchloride/polyvinylchloride_solv_water/checkpoint/polyvinylchloride_solv_water_checkpoint.pkl'))

print(pdir.info, pdir.checkpoint_path)
# load_mol_and_topo(pdir.info.structure_file, pdir.info.monomer_file)
load_mol_and_topo(solvent_pdb, solvent_json)
