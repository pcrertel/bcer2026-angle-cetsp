#!/bin/bash
# build.sh -- compiles the ALNS-BCD into ./alns_cetsp.
# Requirement: gfortran (tested with GCC 12.3 and 13.3). The closed-form block update (max_angle_chain.f90)
# is compiled with -O3 -march=native, and the rest of the code with -O3, as in the experiments of the paper.
set -e
cd "$(dirname "$0")"
D=$(mktemp -d)
cp src/* "$D"/
( cd "$D"
  gfortran -O3 -march=native -c max_angle_chain.f90
  gfortran -O3 instance.f90 blockcoordinatedescent.f90 metaheuristic.f90 main.f90 max_angle_chain.o -o alns_cetsp )
cp "$D"/alns_cetsp .
rm -rf "$D"
echo "built: $(pwd)/alns_cetsp"
