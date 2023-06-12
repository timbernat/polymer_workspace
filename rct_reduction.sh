redux_name='water_soluble_reduced'

python -m setup_collection -strct water_soluble_polymers --output_name wasps_full
python -m polymer_master -src Collections/wasps_full redux -pdb "compatible_pdbs_updated/$redux_name/${redux_name}_structures" -mono "compatible_pdbs_updated/$redux_name/${redux_name}_monomers" -N 200 -f paam_modified peg_modified
python -m setup_collection -strct "$redux_name"  --output_name wasps_reduced