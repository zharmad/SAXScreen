#!/bin/bash

function get_unique_entries() {
    filename=$1 ; column=$2
    awk "{print \$$column }" $filename | sort | uniq
}

function get_files_with_entry() {
    filename=$1 ; column=$2 ; value=$3 ; target=$4
    awk " \$$column == $value {print \$$target} " $filename
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

if [[ "$analysis_comparison_fit_metric" == "V_R" ]] || [[ "$analysis_comparison_fit_metric" == "chi_free" ]] ; then
    fit_args="-Dmax $analysis_comparison_Dmax"
else
    fit_args=""
fi

# Obtain the list of possible titration points from the dictionary.

titration_points_list=$(get_unique_entries $titration_dictionary 2 | tr '\n' ' ')
num_points=$(echo $titration_points_list | wc -w)
titr_first=$(echo $titration_points_list | awk '{print $1}')
ID_list=$(get_files_with_entry $titration_dictionary 2 $titr_first 1)
ID_last=$(echo $ID_list | awk '{print $NF}')

[ ! -e $analysis_output_folder ] && mkdir $analysis_output_folder

echo "= = Note: a total of $num_points distinct titration points found."
for titr in $titration_points_list ; do

    if [ ! -e $analysis_output_folder/fitted_${analysis_comparison_fit_metric}_${titr}_matrix.dat ] ; then
        echo "= = Running fit at titration point $titr ..."
        ID_list=$(get_files_with_entry $titration_dictionary 2 $titr 1)

        tag=""
        for i in $ID_list ; do
            tag="$tag $buffer_subtracted_saxs_folder/${i}_${titr}_${buffer_subtracted_saxs_suffix}.dat"
        done
        ls $tag
        python $script_location/fit-saxs-curves.py \
            -metric $analysis_comparison_fit_metric \
            -mode 2 -qmin $q_min -qmax $q_max $fit_args \
            -o $analysis_output_folder/fitted_${analysis_comparison_fit_metric}_${titr} $tag
    else
        echo "= = Skipping titration point as the output file has been found."
    fi

done

# = = = = =
echo "= = = Copying template gnuplot script and customising for this screen..."
echo "= = = NB: Tnis does not deal with missing titration points."

tag2=""
x=0
for j in $ID_list ; do
    tag2="$tag2\"$j\"\\ $x"
    if [[ "$j" != "$ID_last" ]] ; then
        tag2="$tag2, "
    fi
    x=$((x+1))
done

sed "s/VAR_NUM/$((x-1))/;s/VAR_AXESLABELS/$tag2/;s/VAR_IDLIST/$titration_points_list/;s/VAR_METRIC/$analysis_comparison_fit_metric/;" \
    $script_location/plot-chimatrix-template.gnuplot \
    > $analysis_output_folder/chimatrix-plot.gnuplot



