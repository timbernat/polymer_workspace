python transfer_charges.py       -src wasps_reduced -out wasps_full -pdb water_soluble_polymers
python assign_partial_charges.py -src wasps_full -cp long_chain_chg_params -s unsolv
python vacuum_anneal_conf.py     -src wasps_full -sim vacuum_anneal -r 3