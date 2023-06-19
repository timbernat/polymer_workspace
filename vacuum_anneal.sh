#!/bin/bash

# arg passing
coll_dir=$1 # 'Collections'
full_name=$2 # "wasps_full"
conf_name=$3
num_confs=$4 #4

# arg processing
full_path="$coll_dir/$full_name"
conf_path="$coll_dir/$conf_name"

python -m setup_collection --output_path $conf_path # create empty destination Collection
for ((i=1; i <= $num_confs; i++));
do
    echo $i;
    python -m polymer_master -src $full_path  simulate -sp vacuum_anneal
    python -m polymer_master -src $full_path  transfer_struct -dest $conf_path -caf conf_$i
done

python -m polymer_master -src $conf_path  solvate --solvents WATER_TIP3P
python -m polymer_master -src $conf_path  charge  -cp long_chain_chg_params --solvent --no-charge # temporary measure to ensure conformers also have charges assigned