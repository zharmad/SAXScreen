#!/bin/bash

#function get_unique_entries() {
#    filename=$1 ; column=$2
#    awk "{print \$$column }" $filename | sort | uniq
#}

#function get_files_with_entry() {
#    filename=$1 ; column=$2 ; value=$3 ; target=$4
#    awk " \$$column == $value {print \$$target} " $filename
#}

function filter_files() {
    keys=($2)
    targ=($3)
    output=()
    for i in ${!keys[*]} ; do
        if [[ "$1" == "${keys[i]}" ]] ; then
            output+=("${targ[i]}")
        fi
    done
    echo ${output[*]}
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

if [[ "$analysis_comparison_fit_metric" == "V_R" ]] || [[ "$analysis_comparison_fit_metric" == "chi_free" ]] ; then
    fit_args="-Dmax $analysis_comparison_Dmax"
else
    fit_args=""
fi

# = = = Obtain the list of possible titration points from the dictionary.
# Get number of headers and remove the last two columns as they are special.

# Iterate through the titration dictionary can collect data on all files.
nHead=$(count_headers $titration_dictionary)
nHead=$((nHead-2))
colA=$(search_header_column $titration_dictionary ligandName)
colB=$(search_header_column $titration_dictionary ligandReceptorRatio)

titrationFiles=()
ligandNames=()
prevName=""
ligandRatios=()

while read line
do
    [[ "$line" == "" || "${line:0:1}" == "#" ]] && continue
    set -- $line
    # Note: The dictionary is expected to be of form
    eval inputPrefix=$(form_substring $nHead)
    outputPrefix=$inputPrefix

    # = = = Identify buffer-subtracted surve file(s).
    sourceFile=$(ls $buffer_subtracted_saxs_folder/${inputPrefix}_${buffer_subtracted_saxs_suffix}.dat)
    sourceFile=$(pick_if_multiple NF $sourceFile)
    assert_file $sourceFile

    eval ligName=\$$colA
    eval ligRatio=\$$colB
    if [[ "$prevName" != "$ligName" ]] ; then
        ligandNames+=("$ligName")
        prevName=$ligName
    fi
    ligandRatios+=("$ligRatio")
    if [[ "$sourceFile" == "" ]] ; then
        echo "= = WARNING: $sourceFile is missing, although it contains  corresponding titration entry."
        titrationFiles+=("MISSING_SKIP")
    else
        titrationFiles+=("$sourceFile")
    fi

done < $titration_dictionary

# = = = Now with a complete list, obtain only the unique ligand ratios
IFS=$'\n' uniqueRatios=($(sort -u <<< "${ligandRatios[*]}"))
unset IFS

nPoints=${#uniqueRatios[*]}
nLigands=${#ligandNames[*]}
#titrFirst=$(echo $uniqueRatios | awk '{print $1}')
#ID_list=$(get_files_with_entry $titration_dictionary 2 $titr_first 1)
#ID_last=$(echo $ID_list | awk '{print $NF}')

[ ! -e $analysis_output_folder ] && mkdir $analysis_output_folder
echo "= = NOTE: a total of $nPoints distinct ligand ratios found."
echo "    ... ${ligandRatios[*]}"
echo "= = NOTE: a total of $nLigands ligand titrations detected."
echo "    ... ${ligandNames[*]}"

#echo ${titrationFiles[*]}

# = = = For each ligand ratio, create matrix of comparison values

for titr in ${uniqueRatios[*]} ; do
    outputPrefix=fitted_${analysis_comparison_fit_metric}_${titr}
    outputFile=$analysis_output_folder/${outputPrefix}_matrix.dat
    if [ ! -e $outputFile ] ; then
        echo "= = Running fit across all curves at ligand ratio $titr ..."
        fileList=$(filter_files $titr "${ligandRatios[*]}" "${titrationFiles[*]}")
        echo "    ... $fileList"
        python $script_location/fit-saxs-curves.py \
            -metric $analysis_comparison_fit_metric \
            -mode 2 -qmin $q_min -qmax $q_max $fit_args \
            -o $analysis_output_folder/$outputPrefix \
            $fileList
    else
        echo "= = Skipping titration point as the output file exists: $outputFile"
    fi

done

# = = = = =
echo "= = = Fits complete."
echo "= = = Copying template gnuplot script and customising for this screen..."
echo "= = = NB: This does not deal with missing titration points."

tag2=""
for j in ${!ligandNames[*]} ; do
    tag2="$tag2\"${ligandNames[j]}\"\\ $j"
    if [[ $j -lt $((nLigands-1)) ]] ; then
        tag2="$tag2, "
    fi
done

sed "s/VAR_NUM/$j/;s/VAR_AXESLABELS/$tag2/;s/VAR_IDLIST/${uniqueRatios[*]}/;s/VAR_METRIC/$analysis_comparison_fit_metric/;" \
    $script_location/plot-chimatrix-template.gnuplot \
    > $analysis_output_folder/chimatrix-plot.gnuplot



