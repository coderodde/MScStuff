# Scale font and line width (dpi) by changing the size! It will always display stretched.
set terminal svg size 400,300 enhanced fname 'arial'  fsize 10 butt solid
set output 'out.svg'

# Key means label...
set key inside bottom right
set xlabel 'Size (1000)'
set ylabel 'Time (milliseconds)'
set title 'Strongly connected component algorithms'
plot  "data.txt" using 1:2 title 'Kosaraju' with lines, "data.txt" using 1:3 title 'Path-based' with lines, "data.txt" using 1:4 title 'Tarjan' with lines
