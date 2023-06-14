pdb_dir='compatible_pdbs_updated'
coll_dir='Collections'

full_name='water_soluble_polymers'
redux_name='water_soluble_reduced'

python -m setup_collection -pdb "$pdb_dir/$full_name/${full_name}_structures" -mono "$pdb_dir/$full_name/${full_name}_monomers" --output_path "$coll_dir/$full_name"
python -m polymer_master -src "$coll_dir/$full_name" redux -pdb "$pdb_dir/$redux_name/${redux_name}_structures" -mono "$pdb_dir/$redux_name/${redux_name}_monomers" -N 200 -f paam_modified peg_modified
python -m setup_collection -pdb "$pdb_dir/$redux_name/${redux_name}_structures" -mono "$pdb_dir/$redux_name/${redux_name}_monomers" --output_path "$coll_dir/$redux_name"