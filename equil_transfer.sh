coll_dir=$1
conf_name=$2
equil_name=$3

conf_path="$coll_dir/$conf_name"
equil_path="$coll_dir/$equil_name"

python -m setup_collection --output_path $equil_path # create empty destination Collection
python -m polymer_master -src $conf_path  transfer_struct -dest $equil_path -caf equil --solvent
python -m polymer_master -src $equil_path  charge -cp long_chain_chg_params --no-charge # temporary measure to ensure conformers also have charges assigned