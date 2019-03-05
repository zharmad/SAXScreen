#!/bin/bash

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

# ReceptorName ReceptorConc LigandName LigandReceptorRatio LigandBufferRatio FileLocation
function count_headers() {
    head -n 1 $1 | sed 's/[#@%]//g' | awk '{print NF}'
}

function search_header_column() {
    head -n 1 $1 | sed 's/[#@%]//g' | awk -v targ=$2 \
    '{ for (i=1;i<=NF;i++) { if ( tolower($i) == tolower(targ)) { printf "%i", i ; exit } } print -1 }'
}

get_general_parameters

header_query=ligandReceptorRatio
col=$(search_header_column $titration_dictionary $header_query)
echo $col
if [[ $col -eq -1 ]] ; then
    echo " = = ERROR: the header column search from $titration_dictionary failed to match the query: $header_query"
    exit
else
    echo "    ...the titration dictionary gives column $col as the ligand ratio."
fi

# Get the averaged apo-sample from the titration dictionary
echo "= = Averaging the duplicate non-subtracted apo receptor curves..."
[ -e $average_raw_apo_sample ] && rm $average_raw_apo_sample
fileList=$(awk -v col=$col '$col == 0.0 {print $NF}' $titration_dictionary)
echo "    ...files:" $fileList
python $script_location/average-multicurve.py \
    -o $average_raw_apo_sample -s none -f 1 $fileList
if [ ! -e $average_raw_apo_sample ] ; then
    echo "= = ERROR in generation of averaged raw apo receptor curves!"
    exit 1
fi

# Get average across all sample buffers
echo "= = Averaging the duplicate non-subtracted buffer curves..."
[ -e $average_raw_sample_buffer ] && rm $average_raw_sample_buffer
python $script_location/average-multicurve.py \
    -o $average_raw_sample_buffer -s none -f 1 $(cat $sample_buffer_list)
if [ ! -e $average_raw_sample_buffer ] ; then
    echo "= = ERROR in generation of averaged raw apo receptor curves!"
    exit 1
fi

# Obtain the difference between the averages to serve as an overall reference point for the subtracted apo scattering.
echo "= = Performing subtracting between the two averages..."
[ -e $average_subtracted_sample ] && rm $average_subtracted_sample
python $script_location/manual-subtract.py \
    -s $average_raw_apo_sample -o $average_subtracted_sample $average_raw_sample_buffer
if [ ! -e $average_subtracted_sample ] ; then
    echo "= = ERROR in generation of averaged raw apo receptor curves!"
    exit 1
fi

# = = = Following is where ligand scattering has been measured.
bLigand=$(uncased_string_match $use_ligand_scattering yes)
if [[ "$bLigand" == "True" ]] ; then
    echo "= = Performing operations for ligand buffers to obtain the buffer difference curve.."
    # Get average across all ligand buffers
    python $script_location/average-multicurve.py \
        -o $average_raw_ligand_buffer -s none -f 1 $(cat $ligand_buffer_list)

    # Get the difference scattering between the two buffers for later correction.
    python $script_location/manual-subtract.py \
        -s $average_raw_sample_buffer -o $buffer_difference_file $average_raw_ligand_buffer
else
    echo "= = = NOTE: Ligand scattering has not been set to 'yes', and so ignoring ligand-based averaging."
fi


