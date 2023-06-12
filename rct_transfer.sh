full_coll="Collections/wasps_full"
redux_coll="Collection/wasps_reduced"

python -m polymer_master -src $redux_coll transfer -targ $full_coll
python -m polymer_master -src $full_coll  charge   -cp long_chain_chg_params --no-solvent
python -m polymer_master -src $full_coll  anneal   -sim vacuum_anneal -r 3