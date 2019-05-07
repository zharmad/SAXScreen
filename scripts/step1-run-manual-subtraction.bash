#!/bin/bash

# = = This script runs a manusl subtraction process.
#
#  It takes the raw SAXS file and obtains the synchrotron number from the file name,
#  ..then it attempts to resolve the buffer curves directly before and after the measurement.
#  ..then an autocorrection is carried out when there is a difference between the sample and ligand buffer.

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

ofold=$buffer_subtracted_saxs_folder
[ ! -e $ofold ] && mkdir $ofold

# Iterate through the titration dictionary and conduct the buffer subtraction process
while read line
do
    [[ "$line" == "" || "${line:0:1}" == "#" ]] && continue
    set -- $line
    # Note: The dictionary is expected to be of form
    # Ligand_ID  Ligand:Protein_ratio  ligand:sample_raw_fraction File_location
    # The Synchrotron numbering system varies and will be determined by the file name and type.

    nHead=$(count_headers $titration_dictionary)
    nHead=$((nHead-2))
    eval output_prefix=$(form_substring $nHead)
    colB=$(search_header_column $titration_dictionary ligandBufferRatio)
    eval ligBuffRatio=\$$colB
    colB=$(search_header_column $titration_dictionary fileLocation)
    eval fileLoc=\$$colB

    ofile=$ofold/${output_prefix}_${buffer_subtracted_saxs_suffix}.dat
    if [ -e $ofile ] ; then
        echo "= = NOTE: Output file already exists! Skipping. ( $ofile )"
        shift ; shift
        continue
    fi
    echo "= = Running for ${output_prefix}..."
    if (( $(echo "$ligBuffRatio > 0.0" | bc -l) )) ; then
        echo "= = = Correcting for buffer differences between protein and ligand buffer: $ligBuffRatio "
        addstr="--add_file $buffer_difference_file --add_mult $ligBuffRatio "
    else
        echo "= = = No buffer correction as $ligBuffRatio evaluates to zero."
        addstr=""
    fi

    # Construct file names differently
    case $synchrotron_source in
        Grenoble|ESRF)
            # All files named according to ./1d/*_<SynchID>_ave.dat
            input=$fileLoc
            if [ ! -e $input ] ; then
                echo "= = ERROR: The source file has not been found! ( $input ) Aborting."
                exit 1
            fi
            # Figure out Synchrotron ID from file location. Strip path, prefixes and suffixes
            file_path=${input%/*}
            tmp=${input%_ave.dat} ; s=${tmp##*_} ; y=$(echo $s | sed 's/^0*//')
            x=$(printf "%03d" $((y-1)) )
            z=$(printf "%03d" $((y+1)) )
            buffer1=$(ls ${synchrotron_data_location_sample}/*_${x}_ave.dat) ; buffer2=$(ls ${synchrotron_data_location_sample}/*_${z}_ave.dat)
            assert_files $input $buffer1 $buffer2
            python $script_location/manual-subtract.py \
                -o $ofile -s $input $buffer1 $buffer2 $addstr
            if [ ! -e $ofile ] ; then
                echo "= = ERROR: Output file not created! Something seems to be wrong with the input. Aborting."
                exit 1
            fi
            ;;
        Hamburg|DESY)
            echo "= = Please wait until Hamburg data is correctly coded. Aborting."
            exit 2
            ;;
        Diamond|Chilton|Oxford)
            echo "= = WARNING: Manual subtraction not implemented. Copying over synchrotron-supplied subtracted curves."
            input=$( echo $fileLoc | sed 's,av_dat,sub_dat,;s,_av\.dat,_sub.dat,')
            cp $input $ofile
            ;;
        *)
            echo "= = ERROR: Synchtron source not named correctly in setting file? ( $synchrotron_source ) Unable to figure out how to treat raw data files."
            exit 1
            ;;
    esac
done < $titration_dictionary
