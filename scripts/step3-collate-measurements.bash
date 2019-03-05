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

function search_header_column() {
    head -n 1 $1 | sed 's/[#@%]//g' | awk -v targ=$2 \
    '{ for (i=1;i<=NF;i++) { if ( tolower($i) == tolower(targ)) { printf "%i", i ; exit } } print -1 }'
}

function pick_if_multiple() {
    if [[ $# -gt 2 ]] ; then
        field=$1 ; shift
        echo "= = WARNING: more than one file found according to the settings given for file trawling! : $*" > /dev/stderr
        echo "    ...will use the file #$field." > /dev/stderr
        echo $* | awk "{print \$$field}"
    else
        echo $2
    fi
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
        #file_suffix=_${buffer_subtracted_saxs_suffix}.dat
        file_suffix=*
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

    # Get number of headers and remove the last two columns as they are special.
    nHead=$(count_headers $titration_dictionary)
    nHead=$((nHead-2))
    eval input_prefix=$(form_substring $nHead)
    colA=$(search_header_column $titration_dictionary ligandName)
    colB=$(search_header_column $titration_dictionary ligandReceptorRatio)
    if [[ $colA -eq -1 ]] || [[ $colB -eq -1 ]] ; then
        echo " = = ERROR: the header column search from $titration_dictionary failed to match queries for ligandName and ligandReceptorRatio"
        exit 1
    fi
    eval ligName=\$$colA
    eval ligRatio=\$$colB

    # Change to file search to accept wildcard entries in the file settings.
    if ! file=$(ls $file_prefix/$input_prefix$file_suffix) ; then
        echo "= = WARNING: skipping entry because file has not been found!" > /dev/stderr
        continue
    fi
    # Check if more than one file is found:
    if [[ $(echo $file | wc -w) -gt 1 ]] ; then
        echo "= = WARNING: more than one file found according to the settings given for file trawling! : $file" > /dev/stderr
        echo "    ...will use the last file." > /dev/stderr
        file=$(echo $file | awk '{print $NF}')
    fi
    if [[ "$prevID" != "$ligName" ]] ; then
        if [[ "$prevID" != "" ]] ; then
            echo "&"
        fi
        # Start of new titration
        echo "@s$s legend \"$ligName\""
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
    echo $ligRatio $out
    prevID=$ligName

done < $titration_dictionary
echo "&"
} > $ofile
