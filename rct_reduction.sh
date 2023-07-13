pdb_dir='compatible_pdbs_updated'

coll_dir=$1
full_name=$2
redux_name=$3

python -m setup_collection --output_path "$coll_dir/$full_name"  -pdb "$pdb_dir/$full_name/${full_name}_structures"   -mono "$pdb_dir/$full_name/${full_name}_monomers" 
python -m polymer_master   -src "$coll_dir/$full_name" redux     -pdb "$pdb_dir/$redux_name/${redux_name}_structures" -mono "$pdb_dir/$redux_name/${redux_name}_monomers" -N 200 -f peg_modified paam_modified pnipam_modified
python -m setup_collection --output_path "$coll_dir/$redux_name" -pdb "$pdb_dir/$redux_name/${redux_name}_structures" -mono "$pdb_dir/$redux_name/${redux_name}_monomers" 