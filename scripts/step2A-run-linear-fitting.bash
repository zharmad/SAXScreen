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

assert_file $average_raw_apo_sample $average_raw_sample_buffer $average_raw_ligand_buffer

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
    # Ligand_ID  Ligand:Protein_ratio  ligand:sample_raw_fraction  File_location
    target_file=$4

    # Look up the raw scattering file corresponding to the name of the titration
    ligand_file=$(awk -v sp=$1 '$1 == sp {print $NF}' $ligand_dictionary )
    if [[ "$ligand_file" == "" ]] ; then echo "= = WARNING: ligand file not found in dictionary $ligand_dictionary with pattern $1!" ; exit 1 ; fi
    assert_file $target_file $ligand_file

    echo "= = Running linearity fitting for ${1} ${2}"

    outpref=$ofold/${1}_${2}
    [ -e ${outpref}_model.dat ] && continue
    python $script_location/calculate-linearity.py \
        -a $q_min -b $q_max --doError --ntrials $linear_fitting_error_trials \
        -f $target_file -o $outpref \
        $average_raw_apo_sample $ligand_file $average_raw_sample_buffer $average_raw_ligand_buffer

    chi=$(grep chi ${outpref}_model.dat | awk '{print $(NF-2), $NF}')

done < $titration_dictionary

