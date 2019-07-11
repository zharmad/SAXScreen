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
colFile=$(search_header_column $titration_dictionary fileLocation)
apoList=""
while read line
do
    [[ "$line" == "" || "${line:0:1}" == "#" ]] && continue
    set -- $line
    eval ligRatio=\$$colB  
    if (( $(echo "$ligRatio > 0.0" | bc -l) )) ; then continue ; fi
    eval inputPrefix=$(form_substring $nHead)
    eval sourceFile=\$$colFile
    apoList="$apoList $sourceFile" 
done < $titration_dictionary

$ATSAS_location/primus${ATSAS_suffix} $apoList $average_raw_apo_sample
