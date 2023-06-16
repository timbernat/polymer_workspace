coll_dir='Collections'
full_coll="wasps_full"
# full_coll="water_soluble_polymers"

python -m setup_collection --output_path "$coll_dir/${full_name}_confs" # create empty destination Collection
python -m polymer_master -src "$coll_dir/$full_coll"  simulate -sp vacuum_anneal
python -m polymer_master -src "$coll_dir/$full_coll"  transfer_struct -sp vacuum_anneal

python -m polymer_master -src "$coll_dir/$full_coll"  solvate --solvents WATER_TIP3P
python -m polymer_master -src "$coll_dir/$full_coll"  charge  -cp long_chain_chg_params --no-solvent --no-charge # temporary measure to ensure conformers also have charges assigned