unset key
set cblabel "reduced-{/Symbol c}"
set cbrange [0:5]
set palette defined (0 "black", 1 '#2233bb', 2 "yellow",  3 "red", 4 "pink", 5 "white")
#set palette defined (0 "black", 1 "red", 3 "yellow",  5 "blue")
set size square
set yrange  [VAR_NUM.5:-0.5]
set x2range [-0.5:VAR_NUM.5]
set xrange  [-0.5:VAR_NUM.5]
set xtics format ""

set x2tics ( VAR_AXESLABELS ) \
            rotate by 90 left
set ytics  ( VAR_AXESLABELS )

#set cbtics add ("noise" 1)
set margins at screen 0.095, at screen 1.03, at screen 0.02, at screen 0.89

set term pngcairo enhanced color \
    notransparent \
    linewidth 2.0 \
    font "Courier-Bold,34" \
    size 2100,1900

ratiolist='VAR_IDLIST'
do for [ratio in ratiolist] {
    filename='fitted_VAR_METRIC_'.ratio.'_matrix.dat'
    outfn='./fitted_VAR_METRIC_'.ratio.'_matrix.png'
    #outfn='./chiMatrix-linear-fc-'.ratio.'.pdf'
    #outfn='./chiMatrix-linear-fc-'.ratio.'.svg'
    set output outfn

    labelstr1='smp:lig'
    labelstr2='1.0:'.ratio

    set label labelstr1 at screen 0.94, 0.97 center font "Arial,32"
    set label labelstr2 at screen 0.94, 0.94 center font "Arial-Bold,40"

    plot filename matrix with image
    unset label

    pause -1 "Press enter to continue..."
}
show size
