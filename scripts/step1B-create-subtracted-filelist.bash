#!/bin/bash

# This script just outputs a list of subtracted files that will be processed
# into a simple text file for usage later.

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
source $script_location/header_functions.bash

# Iterate through the titration dictionary
{
while read line
do
    [[ "$line" == "" || "${line:0:1}" == "#" ]] && continue
    set -- $line
   
    # ReceptorName ReceptorConc LigandName LigandRatio LigRecRatio fileLocation
    # Get number of headers and remove the last two columns as they are special.
    nHead=$(count_headers $titration_dictionary)
    nHead=$((nHead-2))
    eval input_prefix=$(form_substring $nHead)
    source_file=$(ls $buffer_subtracted_saxs_folder/${input_prefix}_${buffer_subtracted_saxs_suffix}.dat)
    source_file=$(pick_if_multiple NF $source_file)
    assert_files $source_file

    [[ "$source_file" == "" ]] && continue

    echo $source_file

done < $titration_dictionary
} > $buffer_subtracted_saxs_filelist
