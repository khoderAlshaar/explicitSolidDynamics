set terminal pdfcairo size 16cm,6cm font ",9"
set output "L2_residuals_side_by_side.pdf"

stressFile = "log.solidRoeFoam.analyzed/custom02_addrhoutoresiduals"
momentumFile = "log.solidRoeFoam.analyzed/residuals"

set multiplot layout 1,2

set xlabel "Time [s]"
set ylabel "Stress L2 residual"
set grid
set xrange [0.00001:]
plot stressFile using 1:2 with lines lw 1.5 title "Stress L2 residual"

set xlabel "Time [s]"
set ylabel "Linear momentum L2 residual"
set grid
plot momentumFile using 1:2 with lines lw 1.5 title "Momentum L2 residual"

unset multiplot
