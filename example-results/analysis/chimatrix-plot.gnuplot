unset key
set cblabel "reduced-{/Symbol c}"
set cbrange [0:5]
set palette defined (0 "black", 1 '#2233bb', 2 "yellow",  3 "red", 4 "pink", 5 "white")
#set palette defined (0 "black", 1 "red", 3 "yellow",  5 "blue")
set size square
set yrange  [7.5:-0.5]
set x2range [-0.5:7.5]
set xrange  [-0.5:7.5]
set xtics format ""

set x2tics ( "RNA01" 0, "RNA02" 1, "RNA03" 2, "RNA04" 3, "RNA05" 4, "RNA06" 5, "RNA07" 6, "RNA08" 7, "RNA09" 8, "RNA10" 9, "RNA11" 10, "RNA12" 11 ) \
            rotate by 90 left
set ytics  ( "RNA01" 0, "RNA02" 1, "RNA03" 2, "RNA04" 3, "RNA05" 4, "RNA06" 5, "RNA07" 6, "RNA08" 7, "RNA09" 8, "RNA10" 9, "RNA11" 10, "RNA12" 11 )

#set cbtics add ("noise" 1)
set margins at screen 0.095, at screen 1.03, at screen 0.02, at screen 0.89

set term pngcairo enhanced color \
    notransparent \
    linewidth 2.0 \
    font "Courier-Bold,34" \
    size 2100,1900

ratiolist='0.0 0.1 0.2 0.4 0.6 0.8 0.9 1.0 '
do for [ratio in ratiolist] {
    filename='chimatrix-'.ratio.'-chiMatrix.dat'
    outfn='./chimatrix-'.ratio.'-chiMatrix.png'
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
