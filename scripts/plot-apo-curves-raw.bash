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

apoList=$(grep "0.0 0.0" $titration_dictionary | awk '{print $NF}')

$ATSAS_location/primus${ATSAS_suffix} $apoList $average_raw_apo_sample
