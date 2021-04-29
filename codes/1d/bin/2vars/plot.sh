for infile in fort.1????
do
    gnuplot <<EOF
    set terminal gif
    set output "${infile}.gif"
    set size ratio -1
    set grid
    set xrange[0:2]
    set yrange[-1:1]
    set xlabel "x"
    set ylabel "y"
    plot '${infile}'lt 0 w filledcurve
EOF
done
convert -delay 10 -loop 0 *.gif movie.gif
