#!/bin/bash

coll_dir='Collections'
full_name="wasps_full"
conf_name="${full_name}_confs"
# full_coll="water_soluble_polymers"
full_path="$coll_dir/$full_name"
conf_path="$coll_dir/$conf_name"

num_confs=4

python -m setup_collection --output_path $conf_path # create empty destination Collection
for ((i=1; i <= $num_confs; i++));
do
    echo $i;
    python -m polymer_master -src $full_path  simulate -sp vacuum_anneal
    python -m polymer_master -src $full_path  transfer_struct -dest $conf_path -caf conf_$i
done

python -m polymer_master -src $conf_path  solvate --solvents WATER_TIP3P
python -m polymer_master -src $conf_path  charge  -cp long_chain_chg_params --no-solvent --no-charge # temporary measure to ensure conformers also have charges assigned