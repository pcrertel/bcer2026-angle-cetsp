#!/bin/bash
# run.sh -- compares the cyclic and the queue block selection rules of the BCD (Section 5.1 of the
# paper) on the three sets of 100 random instances of instances/. Builds bcd_compare from the
# sources of the ALNS-BCD (../alns/src) and runs run_bcd.py. One process, one core.
set -e
cd "$(dirname "$0")"
T=$(mktemp -d)
cp ../alns/src/* bcd_compare.f90 "$T"/
( cd "$T"
  g++ -O3 -march=native -std=c++17 -fno-math-errno -c chain_c_api.cpp -o chain_c_api.o
  gfortran -O3 instance.f90 blockcoordinatedescent.f90 bcd_compare.f90 chain_c_api.o -lstdc++ -o bcd_compare )
cp "$T"/bcd_compare .
rm -rf "$T"
python3 run_bcd.py
