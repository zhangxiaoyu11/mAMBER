#!/bin/bash

# First the bonds
plotcmd="plot "
num=0
for i in `ls *bondeq`; do
  title=$(head -n 1 $i)
  title=${title//#}
  title=${title//[[:space:]]}
  plotcmd="$plotcmd '$i' title '$title' with points pt 7 ps 1 lc $num,"
  (( num++ ))
done
plotcmd=${plotcmd%,}

if [[ $num -eq 0 ]]; then exit 0; fi

gnuplot -persist << EOF
set key autotitle columnhead
set title "Bond Equilibrium Lengths in Sampled Conformations"
set xlabel "Equilibrium Length (Angstrom)"
set xrange [0:]
unset ytics
unset ylabel
$plotcmd

EOF
# Do the angles
plotcmd="plot "
num=0
for i in `ls *angleq`; do
  title=$(head -n 1 $i)
  title=${title//#}
  title=${title//[[:space:]]}
  plotcmd="$plotcmd '$i' title '$title' with points pt 7 ps 1 lc $num,"
  (( num++ ))
done
plotcmd=${plotcmd%,}

if [[ $num -eq 0 ]]; then exit 0; fi

gnuplot -persist << EOF
set key autotitle columnhead
set title "Angle Equilibrium Values in Sampled Conformations"
set xlabel "Equilibrium Phase (radians)"
unset ytics
unset ylabel
set xrange [0:3.14]
$plotcmd

EOF

# Now the dihedrals
plotcmd="plot "
num=0
for i in `ls *diheq`; do
  title=$(head -n 1 $i)
  title=${title//#}
  title=${title//[[:space:]]}
  plotcmd="$plotcmd '$i' title '$title' with points pt 7 ps 1 lc $num,"
  (( num++ ))
done
plotcmd=${plotcmd%,}

if [[ $num -eq 0 ]]; then exit 0; fi

gnuplot -persist << EOF
set key autotitle columnhead
set title "Dihedral Equilibrium Values in Sampled Conformations"
set xlabel "Equilibrium Phase (radians)"
set xrange [-3.14:3.14]
unset ytics
unset ylabel
$plotcmd

EOF
