coll_dir='Collections'
full_coll="water_soluble_polymers"
redux_coll="water_soluble_reduced"

python -m polymer_master -src "$coll_dir/$redux_coll" transfer -targ "$coll_dir/$full_coll"
python -m polymer_master -src "$coll_dir/$full_coll"  charge   -cp long_chain_chg_params --no-solvent
python -m polymer_master -src "$coll_dir/$full_coll"  anneal   -sim vacuum_anneal -r 3