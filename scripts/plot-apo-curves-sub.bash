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
source $script_location/header_functions.bash

nHead=$(count_headers $titration_dictionary)
nHead=$((nHead-2))
colB=$(search_header_column $titration_dictionary ligandReceptorRatio)

apoList=""
while read line
do
    [[ "$line" == "" || "${line:0:1}" == "#" ]] && continue
    set -- $line
    eval ligRatio=\$$colB  
    if (( $(echo "$ligRatio > 0.0" | bc -l) )) ; then continue ; fi
    eval inputPrefix=$(form_substring $nHead)
    sourceFile=$(ls $buffer_subtracted_saxs_folder/${inputPrefix}_${buffer_subtracted_saxs_suffix}.dat)
    apoList="$apoList $sourceFile" 
done < $titration_dictionary

$ATSAS_location/primus${ATSAS_suffix} $apoList $average_subtracted_sample
