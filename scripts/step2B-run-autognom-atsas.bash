#!/bin/bash

function extract-Pr() {
    sed -n '/P(R)/,$ p' $1 | awk '{if ( $1 !~ /[0-9]/) { print "#" $0 } else {print $0} }'
}

function extract-Iq() {
    sed -n '/J EXP/,/ Distance / p' $1 | awk '$1 != "" {if ( $1 !~ /[0-9]/) { print "#" $0 } else {print $1, $NF} }'
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

autorg=$ATSAS_location/autorg$ATSAS_suffix
datgnom=$ATSAS_location/datgnom$ATSAS_suffix

assert_files $autorg $datgnom

ofold=$autognom_output_folder
if [ ! -e $ofold ] ; then
    mkdir $ofold
fi

# Iterate through the titration dictionary
while read line
do
    [[ "$line" == "" || "${line:0:1}" == "#" ]] && continue
    set -- $line
   
    # ReceptorName ReceptorConc LigandName LigandRatio LigRecRatio fileLocation
    # Get number of headers and remove the last two columns as they are special.
    nHead=$(count_headers $titration_dictionary)
    nHead=$((nHead-2))
    eval input_prefix=$(form_substring $nHead)
    output_prefix=$input_prefix
    source_file=$(ls $buffer_subtracted_saxs_folder/${input_prefix}_${buffer_subtracted_saxs_suffix}.dat)
    source_file=$(pick_if_multiple NF $source_file)
    assert_files $source_file

    [[ "$source_file" == "" ]] && continue

    echo "= = Processing $source_file ..."
    out_final=$ofold/${output_prefix}_${autognom_output_designation}
    [ -e ${out_final}_Iq.dat ] && continue

    # Cleave off low and high-angle components.
    awk "\$1 > $q_min && \$1 < $q_max { print }" $source_file > temp-inp.dat

    # = = Auto-estimate the Radius of Gyration
    rgdat=$($autorg temp-inp.dat -f ssv)
    # = = Run DATGNOM.
    $datgnom -r ${rgdat%% *} -o $out_final.out temp-inp.dat
    if [ ! -e ${out_final}.out ] ; then
        echo "= = WARNING: DATGNOM has failed to generate $out_final.out!"
        continue
    fi

    # = = Extract P(r) from DATGNOM output
    extract-Pr $out_final.out > ${out_final}_Pr.dat

    # = = Convert P(r) to the smoothed-I(q)
    python $script_location/analyse-distribution.py \
            --int_type Pr --integrate --err --qmax $q_max \
            -f ${out_final}_Pr.dat > ${out_final}_Iq.dat

    rm -f temp-inp.dat

done < $titration_dictionary

