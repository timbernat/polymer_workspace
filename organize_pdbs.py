from shutil import move
from pathlib import Path
import pandas as pd


def inventory_struct_mono(struct_dir : Path, mono_dir : Path) -> pd.DataFrame:
    '''Takes a directory of structure pdb files and a directory of monomer json files and compiles
    an inventory of which molecules are present and whether they have a corresponding set of monomers
    
    Returns as a column-wise Pandas dataframe'''
    data_rows = []
    for pdb_path in struct_dir.glob('*.pdb'):
        mol_name = pdb_path.stem
        json_path = mono_dir/f'{mol_name}.json'

        row_dict = {
            'Species' : mol_name,
            'PDB present' : True,
            'JSON present' : json_path.exists()
        }
        data_rows.append(row_dict)

    return pd.DataFrame(data_rows)

def organize_struct_mono_files(dump_dir : Path, mono_src : Path):
    '''For organizing structure and monomer files in a dump directory '''
    for dir in dump_dir.iterdir():
        if dir.is_dir() and (dir != mono_src): # to encompass the possibility that the json dir is in fact inside the dump folder
            struct_dir = dir/f'{dir.name}_structures'
            mono_dir = dir/f'{dir.name}_monomers'

            struct_dir.mkdir(exist_ok=True)
            mono_dir.mkdir(exist_ok=True)

            for pdb in dir.glob('*.pdb'):
                move(pdb, struct_dir)
                if (json_path := mono_src/f'{pdb.stem}.json').exists():
                    move(json_path, mono_dir)
            
            inventory = inventory_struct_mono(struct_dir, mono_dir)
            inventory.to_csv(dir/f'{dir.name}_inventory.csv')

if __name__ == '__main__':
    PARENT = Path('compatible_pdbs_updated')
    JSON_SRC = PARENT / 'json_files'
    organize_struct_mono_files(PARENT, JSON_SRC)