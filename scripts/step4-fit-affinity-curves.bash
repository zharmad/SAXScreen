#!/bin/bash

function assert_file() {
    for i in $* ; do
        [ ! -e $i ] && echo "= = WARNING: File $i is not found!" > /dev/stderr
    done
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

if [ ! $1 ] ; then
    echo "Usage: ./script <Quantity> - such as Vc Pr Rg I0 chi , etc. "
    exit
fi
quant=$1

input_file=$analysis_output_folder/${quant}_summary.xvg
out_prefix=$analysis_output_folder/${quant}_analysis

python $script_location/fit-affinity-curves.py \
    --targ $input_file \
    -o $out_prefix \
    --conc_rec $sample_concentration \
    --errormode $analysis_fit_errormode \
    --ntrials $analysis_fit_trials
