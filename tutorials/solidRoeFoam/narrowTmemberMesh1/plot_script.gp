set terminal pngcairo size 1200,800
set output 'pressure_vs_x.png'
set title 'Pressure Distribution Over Time'
set xlabel 'x'
set ylabel 'Pxx'
set grid
set key outside
plot 'postProcessing/sample/1.00023e-05/circleIntersection_p.xy' using 1:2 with lines title 'Time 1.00023e-05'
, 'postProcessing/sample/2.00047e-05/circleIntersection_p.xy' using 1:2 with lines title 'Time 2.00047e-05'
, 'postProcessing/sample/2.99964e-05/circleIntersection_p.xy' using 1:2 with lines title 'Time 2.99964e-05'
, 'postProcessing/sample/3.99987e-05/circleIntersection_p.xy' using 1:2 with lines title 'Time 3.99987e-05'
, 'postProcessing/sample/5.00011e-05/circleIntersection_p.xy' using 1:2 with lines title 'Time 5.00011e-05'
, 'postProcessing/sample/6.00034e-05/circleIntersection_p.xy' using 1:2 with lines title 'Time 6.00034e-05'
, 'postProcessing/sample/6.99951e-05/circleIntersection_p.xy' using 1:2 with lines title 'Time 6.99951e-05'
