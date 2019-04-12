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

apoList=$(ls $buffer_subtracted_saxs_folder/*0.0_${buffer_subtracted_saxs_suffix}.dat)

$ATSAS_location/primus${ATSAS_suffix} $apoList $average_subtracted_sample
