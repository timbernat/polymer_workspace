coll_dir=$1
full_name=$2
redux_name=$3

python -m polymer_master -src "$coll_dir/$redux_name" transfer_mono -targ "$coll_dir/$full_name"
python -m polymer_master -src "$coll_dir/$full_name"  charge        -cp long_chain_chg_params --no-solvent