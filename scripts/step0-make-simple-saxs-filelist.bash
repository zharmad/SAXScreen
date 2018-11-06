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
case $synchrotron_source in
    Grenoble)
        # For Grenoble the title/code is contained in the first line.
        for file in $(ls $synchrotron_data_location_sample/*_ave.dat) ; do
            line=$(head -n 1 $file)
            if echo $line | grep $grep_args $search_pattern > /dev/null ; then
                title=$(echo $line | sed 's/^.*\[[0-9]\+\] //')
                echo "$title $file" >> $template_1
            else
                if echo $line | grep -v "\[" > /dev/null ; then
                    echo $file >> $sample_buffer_list
                fi
            fi
        done
        for file in $(ls $synchrotron_data_location_ligand/*_ave.dat) ; do
            line=$(head -n 1 $file)
            if echo $line | grep $grep_args $search_pattern > /dev/null ; then
                title=$(echo $line | sed 's/^.*\[[0-9]\+\] //')
                echo "$title $file" >> $ligand_dictionary
            else
                echo $file >> $ligand_buffer_list
            fi
        done
        ;;
    Hamburg)
        # For Hamburg, the title/code is contained on the second line. Sample description is on the first line.
        for file in $(ls $synchrotron_data_location_sample/data_*.dat) ; do
            code=$(awk 'NR==2 {print $NF ; exit }' $file )
            if echo $code | grep $grep_args $search_pattern > /dev/null ; then
                echo "$code $file" >> $template_1
            else
                echo $file >> $sample_buffer_list
            fi
        done
        for file in $(ls $synchrotron_data_location_ligand/data_*.dat) ; do
            code=$(awk 'NR==2 {print $NF ; exit }' $file )
            if echo $code | grep $grep_args $search_pattern > /dev/null ; then
                echo "$code $file" >> $ligand_dictionary
            else
                echo $file >> $ligand_buffer_list
            fi
        done
        ;;
    *)
        echo "= = ERROR: Synchrotron source unknown! ( $synchrotron_source )"
        exit 1
        ;;
esac
