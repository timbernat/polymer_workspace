python transfer_monomers.py      -src wasps_reduced -targ wasps_full
python assign_partial_charges.py -src wasps_full -cp long_chain_chg_params -s unsolv
python vacuum_anneal_conf.py     -src wasps_full -sim vacuum_anneal -r 3