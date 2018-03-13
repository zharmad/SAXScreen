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
    awk '{ if ($1==0) {print $2, $3; exit } }' $1
}

function collect-Rg() {
    grep "Real space:" $1 | awk '{print $5, $7}'
}

function compute-integral() {
    python $script_location/analyse-distribution.py --integrate -f $1 --int_type x --error
}

function compute-PorodInvariant() {
    python $script_location/analyse-distribution.py --integrate -f $1 --int_type x^2 --error
}

function compute-PorodVolume() {
    datporod $1 | awk '{print $1, 0 }'
}

function compute-ParticleVolume() {
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

function compute-UnitlessKratky() {
    buffer=$(guinier-plot -f $1 -qmax $2 -nofile)
    Rg=$(echo "$buffer" | awk '$1 == "Rg" {print $(NF-2)}')
    I0=$(echo "$buffer" | awk '$1 == "I0" {print $(NF-2)}')
    echo "$buffer" > /dev/stderr
    factor=$(= $Rg^2/$I0)
    echo "= = = $Rg ^2 / $I0 equals $factor " > /dev/stderr
    echo "# Using automated Rgyr computation from gnuplot-4 fit using qmax of $2"
    echo "$buffer" | awk '{print "#" $0}'
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
    echo "= = Usage: ./script <Quantity to collect>"
    echo "    ...examples:
    Vc - Volume of Correlation
    I0 - I(0)
    Rg - Radius of Gyration
    Q  - Porod Invariant
    PV - Porod Volume directly from DATPOROD.
    PV2 - Particle Volume from Porod Invariant
    qIq - integral for Vc
    UKr - Unitless Kratky curve. Requires an additional xmax argument to estimare Rgyr.
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
    I0)  cmdPref=collect-I0  ;;
    Rg|Rgyr) cmdPref=collect-Rg ;;
    Q)   cmdPref=compute-PorodInvariant ;;
    PV)  cmdPref=compute-PorodVolume ;;
    PV2)  cmdPref=compute-ParticleVolume  ;;
    qIq) cmdPref=compute-integral  ;;
    UKr)  cmdPref=compute-UnitlessKratk  ;;
    *)
        echo " = = ERROR: argument not recognised! ( $mode ) " > /dev/stderr
        exit 1
        ;;
esac

$cmdPref $file $args
