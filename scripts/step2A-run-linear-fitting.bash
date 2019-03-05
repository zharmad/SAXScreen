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

# ReceptorName ReceptorConc LigandName LigandReceptorRatio LigandBufferRatio FileLocation
function count_headers() {
    head -n 1 $1 | sed 's/[#@%]//g' | awk '{print NF}'
}

# Given an integer, form the string ${1}_${2}_...${N}
function form_substring()  {
    out='$'{1}
    for i in `seq 2 $1` ; do
        out=${out}_'$'{$i}
    done
    echo $out
}

get_general_parameters

if [[ "$use_ligand_scattering" == "yes" ]] ; then
    assert_file $average_raw_apo_sample $average_raw_sample_buffer $average_raw_ligand_buffer
else
    assert_file $average_raw_apo_sample $average_raw_sample_buffer
fi

# The Blank file is not currently used.
ofold=$linear_fit_output_folder
[ ! -e $ofold ] && mkdir $ofold

chifile=$ofold/$linear_fit_chi_summary_file
rm -f $chifile

# Iterate through the titration dictionary
while read line
do
    [[ "$line" == "" || "${line:0:1}" == "#" ]] && continue
    set -- $line
    # Note: The dictionary is expected to be of form

    nHead=$(count_headers $titration_dictionary)
    nHead=$((nHead-2))
    eval input_prefix=$(form_substring $nHead)
    output_prefix=$input_prefix

    target_file="${@: -1}"
    assert_file $target_file
    if [[ "$use_ligand_scattering" == "yes" ]] ; then
        # Look up the raw ligand scattering file corresponding to the name of the titration
        ligand_file=$(awk -v sp=$1 '$1 == sp {print $NF}' $ligand_dictionary )
        if [[ "$ligand_file" == "" ]] ; then
            echo "= = ERROR: ligand file not found in dictionary $ligand_dictionary with pattern $1! Aborting."
            exit 1
        fi
        bufferLigand=$average_raw_ligand_buffer
    else
        ligand_file=""
        bufferLigand=""
    fi

    echo "= = Running linearity fitting for $output_prefix"

    outpref=$ofold/$output_prefix
    [ -e ${outpref}_model.dat ] && continue
    python $script_location/calculate-linearity.py \
        -a $q_min -b $q_max --doError --ntrials $linear_fitting_error_trials \
        -f $target_file -o $outpref \
        $average_raw_apo_sample $ligand_file $average_raw_sample_buffer $bufferLigand

    chi=$(grep chi ${outpref}_model.dat | awk '{print $(NF-2), $NF}')

done < $titration_dictionary

