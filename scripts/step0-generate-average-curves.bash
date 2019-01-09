#!/bin/bash

function get_general_parameters() {
    local settings=./general-settings.txt
    while read line
    do
        [[ "$line" == "" || "${line:0:1}" == "#" ]] && continue
        set -- $line
        eval $1=$2
    done < $settings
}

get_general_parameters

# Get the averaged apo-sample from the titration dictionary
python $script_location/average-multicurve.py \
    -o $average_raw_apo_sample -s none -f 1 $(awk '$2 == 0.0 {print $NF}' $titration_dictionary)

# Get average across all sample buffers
python $script_location/average-multicurve.py \
    -o $average_raw_sample_buffer -s none -f 1 $(cat $sample_buffer_list)

# Obtain the difference between the averages to serve as an overall reference point for the subtracted apo scattering.
python $script_location/manual-subtract.py \
    -s $average_raw_apo_sample -o $average_subtracted_sample $average_raw_sample_buffer

# = = = Following is where ligand scattering has been measured.
if [[ "$use_ligand_scattering" == "yes" ]] ; then
    # Get average across all ligand buffers
    python $script_location/average-multicurve.py \
        -o $average_raw_ligand_buffer -s none -f 1 $(cat $ligand_buffer_list)

    # Get the difference scattering between the two buffers for later correction.
    python $script_location/manual-subtract.py \
        -s $average_raw_sample_buffer -o $buffer_difference_file $average_raw_ligand_buffer
else
    echo "= = = NOTE: Ligand scattering has not been set to 'yes', and so ignoring ligand-based averaging."
fi


