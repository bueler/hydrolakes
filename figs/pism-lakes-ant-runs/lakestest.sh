#!/bin/bash

#(lakestest.sh) test 'lakes' hydrology on Antarctic geometry and constant
#  basal melt rate of 1 cm a-1; compare
#  http://www2.gi.alaska.edu/snowice/glaciers/iceflow/bueler-igs-fairbanks-june2012.pdf

# run preprocess.sh before this script

NN=4

DOIT=
#DOIT=echo

PREFIX=
PISMGO="mpiexec -n $NN ${PREFIX}pismr"

VERTGRID="-Lz 5000 -Lbz 2000 -Mz 51 -Mbz 21 -z_spacing equal"

OPTIONS="-sia_e 5.6 -atmosphere given -atmosphere_given_file pism_Antarctica_5km.nc -surface simple -ocean pik -meltfactor_pik 1.5e-2 -ocean_kill"

HYDRO="-hydrology lakes -hydrology_use_const_bmelt -hydrology_const_bmelt 3.1689e-10 -hydrology_hydraulic_conductivity 1.0e-3"

DURATION=20000

dorun () {
  GRID=$1
  LABEL=$2

  #(from antspinCC.sh)  bootstrapping plus short SIA run for 100 years
  cmd="$PISMGO -skip -skip_max 10 -boot_file pism_Antarctica_5km.nc $GRID $VERTGRID $OPTIONS -y 100 -o pre${LABEL}.nc"
  $DOIT $cmd

  #hydrology only for $DURATION years
  cmd="$PISMGO -i pre${LABEL}.nc $OPTIONS $HYDRO -no_mass -no_energy -no_sia -max_dt 10.0 -y $DURATION -o lakes${LABEL}.nc"
  $DOIT $cmd
}

HUNDREDKMGRID="-Mx 60 -My 60"
FIFTYKMGRID="-Mx 120 -My 120"
TWENTYFIVEKMGRID="-Mx 240 -My 240"
FIFTEENKMGRID="-Mx 400 -My 400"
TENKMGRID="-Mx 600 -My 600"
FIVEKMGRID="-Mx 1200 -My 1200"

# first three regenerate results from IGS talk:
#dorun "$HUNDREDKMGRID" 100km
#dorun "$FIFTYKMGRID" 50km
#dorun "$TWENTYFIVEKMGRID" 25km

# these are more expensive, naturally
dorun "$FIFTEENKMGRID" 15km
#dorun "$TENKMGRID" 10km
#dorun "$FIVEKMGRID" 5km

