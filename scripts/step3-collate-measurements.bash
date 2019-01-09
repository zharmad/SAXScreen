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
    echo "Usage: ./script <Quantity> - such as Vc Pr Rg I0 chi PV PV2 , etc. "
    exit
fi
quant=$1

bError=True
case $quant in
    chi)
        ofile=$analysis_output_folder/chi_summary.xvg
        file_prefix=$linear_fit_output_folder
        file_suffix=_chis.dat
        ;;
    I0)
        ofile=$analysis_output_folder/I0_summary.xvg
        file_prefix=$autognom_output_folder
        file_suffix=_${autognom_output_designation}_Iq.dat
        ;;
    Rg|Rgyr)
        ofile=$analysis_output_folder/Rg_summary.xvg
        file_prefix=$autognom_output_folder
        file_suffix=_${autognom_output_designation}.out
        ;;
    Vc)
        ofile=$analysis_output_folder/Vc_summary.xvg
        file_prefix=$autognom_output_folder
        file_suffix=_${autognom_output_designation}_Iq.dat
        ;;
    VR)
        ofile=$analysis_output_folder/VR_summary.xvg
        file_prefix=$buffer_subtracted_saxs_folder
        file_suffix=_${buffer_subtracted_saxs_suffix}.dat
        bError=False
        ;;
    PV)
        ofile=$analysis_output_folder/PV_summary.xvg
        file_prefix=$autognom_output_folder
        file_suffix=_${autognom_output_designation}.out
        ;;
    PV2)
        ofile=$analysis_output_folder/PV2_summary.xvg
        file_prefix=$autognom_output_folder
        file_suffix=_${autognom_output_designation}_Iq.dat
        ;;
    *)
        echo " = = ERROR: argument not recognised! ( $quant ) " > /dev/stderr
        exit 1
        ;;
esac

ofold=$analysis_output_folder
if [ ! -e $ofold ] ; then
    mkdir $ofold
fi

s=0
prevID=""
# Iterate through the titration dictionary again, with subtle variations to detect new titrations.
{
while read line
do
    [[ "$line" == "" || "${line:0:1}" == "#" ]] && continue
    set -- $line
    # Note: The dictionary is expected to be of form
    # Ligand_ID  Ligand:Protein_ratio  ligand:sample_raw_fraction File_location

    file=$file_prefix/${1}_${2}$file_suffix
    if [[ "$prevID" != "$1" ]] ; then
        if [[ "$prevID" != "" ]] ; then
            echo "&"
        fi
        # Start of new titration
        echo "@s$s legend \"$1\""
        if [[ "$bError" == "True" ]] ; then
            echo "@type xydy"
        else
            echo "@type xy"
        fi
        s=$((s+1))
    fi
    if [[ "$quant" == "VR" ]] ; then
        args=$average_subtracted_sample
    else
        args=""
    fi
    echo "= = Executing: $script_location/compute-quantity.bash $file $quant $args" > /dev/stderr
    out=$($script_location/compute-quantity.bash $file $quant $args)
    echo $2 $out
    prevID=$1
done < $titration_dictionary
echo "&"
} > $ofile
