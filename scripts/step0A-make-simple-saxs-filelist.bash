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

function uncased_string_match() {
    awk -v a=$1 -v b=$2 'BEGIN {
    if ( tolower(a) == tolower(b) ) { print "True" } else { print "False" }
    }'
}

# = = = Most synchrotron use a numbering system as follows:
# <title>_####.dat
# This is leveraged to grab the unique number.
function grab_scattering_number() {
    tmp=${1%.dat}
    echo ${tmp##*_} | sed 's/^[0]\+//'
}

get_general_parameters

if [ ! $1 ] ; then
    echo "= = Usage: ./script [-not] <Title of Sample/Ligands synchrotron> [More Titles]"
    echo "= = This script searches the first lines of the raw measurement files for the Name, and copies ones that satisfy the criteria into lists"
    echo "= = Its aim is to help users to set up the dictionaries. Each argument is parsed as an OR to the grep command."
    echo "= = If the title is found, the file is listed as a sample. At Grenoble, if there is no '[' character it is listed as a buffer."
    exit
fi

if [[ "$1" == "-not" ]] ; then
    grep_args="-v"
    shift
fi
search_pattern=$(echo $* | sed 's/ /\\|/g')
template_1=${titration_dictionary}_template
rm -f $template_1 $sample_buffer_list $ligand_buffer_list $ligand_dictionary

echo "# ReceptorName ReceptorConc LigandName LigandReceptorRatio LigandBufferRatio FileLocation" >> $template_1

echo "= = = Given data location:  $synchrotron_data_location_sample"
sleep 1.5

function debug_fileCount() {
    echo "    ...total files counted: $fileCount"
    echo "    ...sample files counted: $sampleCount"
    echo "    ...buffer files counted: $bufferCount"
    echo "    ...ignored files counted: $ignoreCount"
}
fileCount=0 ; sampleCount=0 ; bufferCount=0 ; ignoreCount=0
case $synchrotron_source in
    Grenoble)
        # For Grenoble the title/code is contained in the first line.
        echo "= = Conducting sample file list generation..."
        for file in $(ls $synchrotron_data_location_sample/*_ave.dat) ; do
            filCount=$((fileCount+1))
            line=$(head -n 1 $file)
            if echo $line | grep $grep_args $search_pattern > /dev/null ; then
                title=$(echo $line | sed 's/^.*\[[0-9]\+\] //')
                echo "$title 100 lig 1.0 0.0 $file" >> $template_1
                sampleCount=$((sampleCount+1))
            elif echo $line | grep -v "\[" > /dev/null ; then
                echo $file >> $sample_buffer_list
                bufferCount=$((bufferCount+1))
            else
                echo "    ...NB: file $file has been ignored."
                ignoreCount=$((ignoreCount+1))
            fi
        done
        debug_fileCount

        bLigand=$(uncased_string_match $use_ligand_scattering yes)
        if [[ "$bLigand" == "True" ]] ; then
          echo "= = Conducting ligand file list generation..."
          for file in $(ls $synchrotron_data_location_ligand/*_ave.dat) ; do
            fileCount=$((fileCount+1))
            line=$(head -n 1 $file)
            if echo $line | grep $grep_args $search_pattern > /dev/null ; then
                sampleCount=$((sampleCount+1))
                title=$(echo $line | sed 's/^.*\[[0-9]\+\] //')
                echo "$title $file" >> $ligand_dictionary
            else
                bufferCount=$((bufferCount+1))
                echo $file >> $ligand_buffer_list
            fi
          done
          debug_fileCount
        else
            echo "= = Skipping ligand file list generation."
        fi

        ;;
    Hamburg|DESY)
        # For Hamburg, the title/code is contained on the second line. Sample description is on the first line.
        # The quirk here is that both sample and buffer files will contain the string used to label buffer measurements. Therefore, test whether it is a sample first.
        echo "= = Conducting sample file list generation..."
        for file in $(ls $synchrotron_data_location_sample/data_*.dat) ; do
            fileCount=$((fileCount+1))
            code=$(awk 'NR==2 {print $NF ; exit }' $file )
            if echo $code | grep $grep_args $search_pattern > /dev/null ; then
                sampleCount=$((sampleCount+1))
                echo "$code 100 lig 1.0 0.0 $file" >> $template_1
            elif echo $code | grep $grep_args buffer > /dev/null ; then
                bufferCount=$((bufferCount+1))
                echo $file >> $sample_buffer_list
            else 
                echo "    ...ignoring file $file , failed to match search pattern"
                ignoreCount=$((ignoreCount+1))
            fi
        done
        debug_fileCount

        bLigand=$(uncased_string_match $use_ligand_scattering yes)
        if [[ "$bLigand" == "True" ]] ; then
          echo "= = Conducting ligand file list generation..."
          for file in $(ls $synchrotron_data_location_ligand/data_*.dat) ; do
            fileCount=$((fileCount+1))
            code=$(awk 'NR==2 {print $NF ; exit }' $file )
            if echo $code | grep $grep_args $search_pattern > /dev/null ; then
                sampleCount=$((sampleCount+1))
                echo "$code $file" >> $ligand_dictionary
            else
                bufferCount=$((bufferCount+1))
                echo $file >> $ligand_buffer_list
            fi
          done
          debug_fileCount
        else
            echo "= = Skipping ligand file list generation."
        fi
        ;;
    Melbourne|Australia)
        # Unique to Australia:the automated averages are somtimes duplicated due to a quirk in the automated workflow. Therefore, check if there is a duplicate. Only store the last one.
        # NB: This will sometimes remove two measurements that have been named identically. Therefore, also check if the numbers are sepated by less than 10. Australia's numbering system is according to raw frames only. 99% of the time independent measurements will be at least 30 apart.
        echo "= = Conducting sample file list generation..."
        echo " = = = Australian Synchrotron standards: "
        echo "    ...using averaged frames in avg/ "
        echo "    ...the first line contains the filename and included averaging frames (not used here)."
        echo "    ...treating all files with name *buffer* as sample buffer files."
        fileRoot=${synchrotron_data_location_sample}/avg
        for file in $(ls $fileRoot/*.dat) ; do
            fileCount=$((fileCount+1))
            currNum=$(grab_scattering_number $file)
            dupSearch=$(ls ${file%????.dat}????.dat)
            if [[ $(echo $dupSearch | wc -w) -gt 1 ]] ; then
                last=$(echo $dupSearch | awk '{print $NF}')
                lastNum=$(grab_scattering_number $last)
                sep=$(awk "BEGIN {print $lastNum - $currNum}" | sed 's/-//')
                #echo $last $file
                if [[ "$last" != "$file" ]] && [[ $sep -lt 10 ]] ; then
                    # This one is a duplicate. Ignore
                    ignoreCount=$((ignoreCount+1))
                    echo "    ...ignoring file $file , duplicate of existing measurement."
                    continue
                fi
            fi
            #echo "= = Processing file: $file"

            # Check if buffer. Get root file name
            header=$(head -n 1 $file)
            if echo $header | grep $grep_args $search_pattern > /dev/null ; then
                sampleCount=$((sampleCount+1))
                title=${header%%_*}
                echo "$title 100 lig 1.0 0.0 $file" >> $template_1
            elif [[ $header == *"buffer"* ]] ; then
                bufferCount=$((bufferCount+1))
                echo $file >> $sample_buffer_list
            else
                ignoreCount=$((ignoreCount+1))
                echo "    ...ignoring file $file , failed to match search pattern."
            fi
        done
        debug_fileCount
        echo "= = NB: Ligand scattering is not yet coded for this source."
        ;;

    *)
        echo "= = ERROR: Synchrotron source unknown! ( $synchrotron_source )"
        exit 1
        ;;
esac
