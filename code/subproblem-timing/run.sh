#!/bin/bash
# run.sh -- time of one closed-form block update (Section 5.1 of the paper), over one million random
# configurations. Builds subproblem_timing with the block update of the ALNS-BCD (../alns/src) and runs it.
# Run it on an idle machine; taskset fixes the core.
set -e
cd "$(dirname "$0")"
T=$(mktemp -d)
cp ../alns/src/max_angle_chain.f90 subproblem_timing.f90 "$T"/
( cd "$T"
  gfortran -O3 -march=native -c max_angle_chain.f90
  gfortran -O3 -march=native subproblem_timing.f90 max_angle_chain.o -o subproblem_timing )
cp "$T"/subproblem_timing .
rm -rf "$T"
taskset -c 0 ./subproblem_timing 1000000 20260901
