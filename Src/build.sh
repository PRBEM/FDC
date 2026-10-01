#!/bin/bash
# usage: ./build.sh <name> [extra flags...]
# e.g.   ./build.sh oldalgo -DROOTS_POLY -DINTN_ORIG -DNO_THETA_SKIP -DOLD_BM
set -e
name=$1; shift
mpif90 -cpp -O2 "$@" *.f90 -I../whamp/src -o "main_$name" \
    -I/usr/lib/x86_64-linux-gnu/openmpi/include -ffree-line-length-none \
    -L/usr/lib/ -l:libblas.so.3 -l:liblapack.so.3 ../whamp/src/libwhamp.a
echo "main_$name <- $(git rev-parse --short HEAD) $(git diff --quiet || echo '+local changes') flags: $*"
