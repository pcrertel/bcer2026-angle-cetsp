#!/bin/bash
# build.sh -- compiles the ALNS-BCD into ./alns_cetsp.
# Requirements: g++ with C++17 and gfortran (tested with GCC 12.3 and 13.3).
#
#   bash build.sh           default build, with std::atan2 (the build used in the paper)
#   bash build.sh --sleef   faster build, with the atan2 of SLEEF (requires libsleef; the results
#                           agree to rounding but are not bit-identical to the default build)
set -e
cd "$(dirname "$0")"
CXXFLAGS="-O3 -march=native -std=c++17 -fno-math-errno"
LIBS=""
if [ "$1" = "--sleef" ]; then
  CXXFLAGS="$CXXFLAGS -DMAXANGLE_USE_SLEEF -include sleef_compat.hpp"
  LIBS="-lsleef"
fi
D=$(mktemp -d)
cp src/* "$D"/
( cd "$D"
  g++ $CXXFLAGS -c chain_c_api.cpp -o chain_c_api.o
  gfortran -O3 instance.f90 blockcoordinatedescent.f90 metaheuristic.f90 main.f90 chain_c_api.o \
     -lstdc++ $LIBS -o alns_cetsp )
cp "$D"/alns_cetsp .
rm -rf "$D"
echo "built: $(pwd)/alns_cetsp"
