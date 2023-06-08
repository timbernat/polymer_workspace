python setup_collection.py       -strct water_soluble_polymers --output_name wasps_full
python reduce_structs_linear.py  -src   wasps_full           --struct_output water_soluble_reduced -N 200 -f paam_modified peg_modified
python setup_collection.py       -strct water_soluble_reduced  --output_name wasps_reduced