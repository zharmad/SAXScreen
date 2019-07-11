#!/bin/bash

function determine_xmgrace_xydy() {
    if grep "xydy" $1 > /dev/null ; then
        echo True
    else
        echo False
    fi
}

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

if [ ! $1 ] ; then
    echo "Usage: ./script <Quantity> - such as Vc Pr Rg I0 chi , etc. "
    exit
fi
quant=$1

input_file=$analysis_output_folder/${quant}_summary.xvg
out_prefix=$analysis_output_folder/${quant}_analysis

# = = = Determine if the file has errorbars.
bError=$(determine_xmgrace_xydy $input_file)
[ $analysis_fit_reject_threshold ] && args="--rejectThreshold $analysis_fit_reject_threshold" \
    || { args="" ; echo "= = NOTE: the rejection threshold has not been set, using script defaults."; }
if [[ "$bError" == "False" ]] || [[ "$analysis_fit_errormode" == "1-point" ]] ; then
    python $script_location/fit-affinity-curves.py \
        --targ $input_file \
        -o $out_prefix \
        --conc_rec $sample_concentration \
        --errormode 1-point $args
else
    [ $analysis_fit_reject_threshold ] && args="--rejectThreshold $analysis_fit_reject_threshold" || { args="" ; echo "= = NOTE: the rejection threshold has not been set, using script defaults."; }
    python $script_location/fit-affinity-curves.py \
        --targ $input_file \
        -o $out_prefix \
        --conc_rec $sample_concentration \
        --errormode $analysis_fit_errormode \
        --nTrials $analysis_fit_trials $args
fi
