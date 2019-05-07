#!/bin/bash

function assert_files() {
    for i in $* ; do
        [ ! -e $i ] && echo "= = WARNING: File $i is not found!" > /dev/stderr
    done
}

# = = = Most synchrotron use a numbering system as follows:
# <title>_####.dat
# This is leveraged to grab the unique number.
function grab_scattering_number() {
    tmp=${1%.dat}
    echo ${tmp##*_} | sed 's/^[0]\+//'
}

# ReceptorName ReceptorConc LigandName LigandReceptorRatio LigandBufferRatio FileLocation
function count_headers() {
    head -n 1 $1 | sed 's/[#@%]//g' | awk '{print NF}'
}

# Given an integer, form the string ${1}_${2}_...${N}
function form_substring()  {
    out='$'{1}
    for i in `seq 2 $1` ; do out=${out}_'$'{$i} ; done
    echo $out
}

# Processes the arguments with unknown number of fields,
# in order to guarantee only a single argument remains.
# If there are more than a single field, then select field 1.
function pick_if_multiple() {
    if [[ $# -gt 2 ]] ; then
        field=$1 ; shift
        echo "= = WARNING: more than one file found in a wildcard file search query : $*" > /dev/stderr
        if [[ "$field" == "NF" ]] ; then
            echo "    ...will use the last file in list." > /dev/stderr
        else
            echo "    ...will use the file #$field." > /dev/stderr
        fi
        echo $* | awk "{print \$$field}"
    else
        echo $2
    fi
}





# Search the first line of data file looking for string matches.
# Then retrun the relevant column entry 
function search_header_column() {
    head -n 1 $1 | sed 's/[#@%]//g' | awk -v targ=$2 \
    '{ for (i=1;i<=NF;i++) { if ( tolower($i) == tolower(targ)) { printf "%i", i ; exit } } print -1 }'
}

function uncased_string_match() {
    awk -v a=$1 -v b=$2 'BEGIN {
    if ( tolower(a) == tolower(b) ) { print "True" } else { print "False" }
    }'
}
