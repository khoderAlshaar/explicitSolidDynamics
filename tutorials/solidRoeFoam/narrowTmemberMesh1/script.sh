#!/bin/bash

# Directory containing the postProcessing results
DATA_DIR="postProcessing/sample"

# Check if the listTimes file exists
if [[ ! -f listTimes ]]; then
    echo "Error: listTimes file not found!"
    exit 1
fi

# Create a temporary gnuplot script
PLOT_SCRIPT="plot_script.gp"
echo "set terminal pngcairo size 1200,800" > $PLOT_SCRIPT
echo "set output 'pressure_vs_x.png'" >> $PLOT_SCRIPT
echo "set title 'Pressure Distribution Over Time'" >> $PLOT_SCRIPT
echo "set xlabel 'x'" >> $PLOT_SCRIPT
echo "set ylabel 'Pxx'" >> $PLOT_SCRIPT
echo "set grid" >> $PLOT_SCRIPT
echo "set key outside" >> $PLOT_SCRIPT
echo -n "plot " >> $PLOT_SCRIPT

# Read the time directories from listTimes and append to the gnuplot script
first=true
while read -r time; do
    FILE="$DATA_DIR/$time/circleIntersection_p.xy"
    if [[ -f $FILE ]]; then
        if $first; then
            echo "'$FILE' using 1:2 with lines title 'Time $time'" >> $PLOT_SCRIPT
            first=false
        else
            echo ", '$FILE' using 1:2 with lines title 'Time $time'" >> $PLOT_SCRIPT
        fi
    else
        echo "Warning: File $FILE not found, skipping..."
    fi
done < listTimes

# Run gnuplot
gnuplot $PLOT_SCRIPT

# Clean up
echo "Plot saved as pressure_vs_x.png"
