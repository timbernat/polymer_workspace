python setup_collection.py       -src water_soluble_polymers --coll_name wasps_full
python create_reduced_structs.py -src wasps_full --output_name water_soluble_reduced -N 200
python setup_collection.py       -src water_soluble_reduced --coll_name wasps_reduced