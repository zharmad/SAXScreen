#!/bin/bash

function extract-Pr() {
    sed -n '/P(R)/,$ p' $1 | awk '{if ( $1 !~ /[0-9]/) { print "#" $0 } else {print $0} }'
}

function extract-Iq() {
    sed -n '/J EXP/,/ Distance / p' $1 | awk '$1 != "" {if ( $1 !~ /[0-9]/) { print "#" $0 } else {print $1, $NF} }'
}

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

autorg=$ATSAS_location/autorg$ATSAS_suffix
datgnom=$ATSAS_location/datgnom$ATSAS_suffix

ofold=$autognom_output_folder
if [ ! -e $ofold ] ; then
    mkdir $ofold
fi

# Iterate through the titration dictionary
while read line
do
    [[ "$line" == "" || "${line:0:1}" == "#" ]] && continue
    set -- $line
    # Note: The dictionary is expected to be of form
    # Ligand_ID  Ligand:Protein_ratio  ligand:sample_raw_fraction File_location

    source_file=$(ls $buffer_subtracted_saxs_folder/${1}_${2}_man01.dat )
    assert_file $source_file
    [[ "$source_file" == "" ]] && continue

    echo "= = Processing $source_file ..."
    out_final=$ofold/${1}_${2}_${autognom_output_designation}

    # Cleave off low and high-angle components.
    awk "\$1 > $q_min && \$1 < $q_max { print }" $source_file > temp-inp.dat

    # Auto-estimate the Radius of Gyration
    rgdat=$($autorg temp-inp.dat -f ssv)
    # Run DATGNOM.
    $datgnom -r ${rgdat%% *} -o $out_final.out temp-inp.dat

    # Extract P(r) from DATGNOM output
    extract-Pr $out_final.out > ${out_final}_Pr.dat

    # Convert P(r) to the smoothed-I(q)
    python $script_location/analyse-distribution.py \
            --int_type Pr --integrate --err \
            -f ${out_final}_Pr.dat > ${out_final}_Iq.dat

    rm -f temp-inp.dat

done < $titration_dictionary
