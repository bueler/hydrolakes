#!/bin/bash

# run newtonconverge.m in matlab/octave after this one

NN=4
CMD="mpiexec -n $NN ./porous -por_steps 1 -por_converge_check"

touch converge.txt
rm converge.txt
for (( TE=40; TE<160; TE=TE+5 )) ; do
    $CMD -da_grid_x 21 -da_grid_y 21 -por_tend_days $TE >> converge.txt
done
for (( TE=12; TE<34; TE=TE+2 )) ; do
    $CMD -da_grid_x 81 -da_grid_y 81 -por_tend_days $TE >> converge.txt
done
for TE in 10.5 11.0 11.5 12.0 12.5 13.0 13.5 14.0 14.5 15.0 ; do
    $CMD -da_grid_x 321 -da_grid_y 321 -por_tend_days $TE >> converge.txt
done

