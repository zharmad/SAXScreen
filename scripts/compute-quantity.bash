#!/bin/bash

# Jochen Hub's command line calculator (use square brackts to avoid bash confusion)
=() {
    local in="$(echo "$@" | sed -e 's/\[/(/g' -e 's/\]/)/g')"
    # echo "in=$in"
    awk 'BEGIN {print '"$in"'}' < /dev/null
}

function collect-chi() {
    awk 'NR == 1 { print $(NF-2), $NF ; exit }' $1
}

function collect-I0() {
    # = = Default use: get value from smoothed I(q) curve. Otherwise call autorg to evalute curve.
    out=$(awk '{ if ($1==0.0) {print $2, $3; exit } }' $1)
    [[ "$out" != "" ]] && echo $out || autorg $1 | awk '$1 == "I(0)" {print $3, $5 ; exit}'
}

function collect-Rg() {
    # = = Default use: get value from GNOM output. Otherwise call autorg to evaluate curve.
    out=$(grep "Real space Rg:" $1 | awk '{print $(NF-2), $NF}')
    [[ "$out" != "" ]] && echo $out || autorg $1 | awk '$1 == "Rg" {print $3, $5 ; exit}'
}

function compute-integral() {
    python $script_location/analyse-distribution.py --integrate -f $1 --int_type x --error
}

function compute-Rg() {
#    autorg -o temp.file $1
#    rm -f temp.file
    autorg $1 | grep Rg | awk '{print $(NF-3), $(NF-1)}'
}

function compute-PorodInvariant() {
    python $script_location/analyse-distribution.py --integrate -f $1 --int_type x^2 --error
}

function compute-PorodVolume() {
    # = = default use on GNOM output. Otherwise call autorg to provide estimates.
    if [[ "${1##*.}" == "out" ]] ; then
        datporod $1 | awk '{print $(NF-1)}'
    else
        read -r Rg I0 <<<$(autorg $1 | awk 'NR==1 {Rg=$3} NR==2 {I0=$3} END {print Rg, I0}')
        datporod --rg=$Rg --i0=$I0 $1 | awk '{print $(NF-1)}'
    fi
}

function compute-ParticleVolume() {
    # = = Do not use. Is unreliable.
    factor=$(= 2.0*3.1415927^2)
    tmp1=$(collect-I0 $1)
    tmp2=$(compute-PorodInvariant $1)
    I0=${tmp1% *}  ; dI0=${tmp1#* }
    Int=${tmp2% *} ; dInt=${tmp2#* }
    Val=$(= $factor*$I0^2/$Int)
    dVal=$(= sqrt[[$dI0/$I0]^2+[$dInt/$Int]^2]*$Val)
    echo $Val $dVal
}

function compute-Vc() {
    tmp1=$(collect-I0 $1)
    tmp2=$(compute-integral $1)
    I0=${tmp1% *}  ; dI0=${tmp1#* }
    Int=${tmp2% *} ; dInt=${tmp2#* }
    Val=$(= $I0/$Int)
    dVal=$(= sqrt[[$dI0/$I0]^2+[$dInt/$Int]^2]*$Val)
    echo $Val $dVal
}

function compute-VR() {
    [ -e temp-Exp-1.xvg ] && rm temp-Exp-1.xvg
    python $script_location/fit-saxs-curves.py \
        -metric V_R -Dmax $analysis_comparison_Dmax \
        -o temp \
        -qmin $q_min -qmax $q_max \
        $1 $2 > /dev/stderr
    if [ -e temp-Exp-1.xvg ] ; then
        grep volatRat temp-Exp-1.xvg | awk '{print $NF}'
        rm -f temp-Exp-1.xvg
    else
        echo "= = ERROR in compute-VR from compute-quantity.bash! The input to the function does not generate an output file!" > /dev/stderr
    fi
}

function compute-UnitlessKratky() {
    #buffer=$(guinier-plot -f $1 -qmax $2 -nofile)
    #Rg=$(echo "$buffer" | awk '$1 == "Rg" {print $(NF-2)}')
    #I0=$(echo "$buffer" | awk '$1 == "I0" {print $(NF-2)}')
    #echo "$buffer" > /dev/stderr
    #factor=$(= $Rg^2/$I0)
    #echo "= = = $Rg ^2 / $I0 equals $factor " > /dev/stderr
    #echo "# Using automated Rgyr computation from gnuplot-4 fit using qmax of $2"
    #echo "$buffer" | awk '{print "#" $0}'
    read -r Rg I0 factor <<<$(autorg $1 | awk 'NR==1 {Rg=$3} NR==2 {I0=$3} END {print Rg, I0, Rg*Rg/I0}')
    python $script_location/analyse-distribution.py --xscale $Rg --yscale $factor --transform -f $1 --int_type x^2 --error
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

if [ ! $2 ] ; then
    echo "= = Usage: ./script <Curve> <Quantity to collect> [Additional Arguments]"
    echo "    ...examples:
    Vc - Volume of Correlation.
    VR - Volatility of Ratio (minimized). Fits VR ( I_A+c , I_B-c ) over constant scattering c to remove potential buffer differences.
    I0 - I(0)
    Rg - Radius of Gyration from GNOM.
    autorg - Radius of Gyration from ATSAS autorg utility.
    Q  - Porod Invariant. Computes just the q^2I(q) integral without considering complications.
    PV - Porod Volume directly from DATPOROD. Accepts either P(r) via GNOM output with *.out extension, or I(q) as anything else.
    PV2 - Particle Volume from Porod Invariant, without considering complications.
    qIq - integral for Vc
    UKr - Unitless Kratky curve. Requires an additional xmax argument to estimate Rgyr.
    "
    exit
fi

file=$1
mode=$2
shift ; shift
args=$*

case $mode in
    chi) cmdPref=collect-chi ;;
    Vc)  cmdPref=compute-Vc  ;;
    VR)  cmdPref=compute-VR  ;;
    I0)  cmdPref=collect-I0  ;;
    Rg|Rgyr) cmdPref=collect-Rg ;;
    autorg) cmdPref=compute-Rg ;;
    Q)   cmdPref=compute-PorodInvariant ;;
    PV)  cmdPref=compute-PorodVolume ;;
    PV2)  cmdPref=compute-ParticleVolume  ;;
    qIq) cmdPref=compute-integral  ;;
    UKr)  cmdPref=compute-UnitlessKratky  ;;
    *)
        echo " = = ERROR: argument not recognised! ( $mode ) " > /dev/stderr
        exit 1
        ;;
esac
$cmdPref $file $args

