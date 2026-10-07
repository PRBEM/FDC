#!/bin/bash
# usage: ./build.sh <name> [extra flags...]
# e.g.   ./build.sh oldalgo -DROOTS_POLY -DINTN_ORIG -DNO_THETA_SKIP -DOLD_BM
#        ./build.sh prof -DPROFILE        ./build.sh diag -DDIAG
set -e
name=$1; shift
FLAGS="-cpp -O2 $* -ffree-line-length-none"
LIBS="-L/usr/lib/ -l:libblas.so.3 -l:liblapack.so.3 ../whamp/src/libwhamp.a"
INC="-I../whamp/src -I/usr/lib/x86_64-linux-gnu/openmpi/include"

# whamp is a submodule; FDC links libwhamp.a and uses its .mod files
if [ ! -f ../whamp/src/libwhamp.a ] || [ ! -f ../whamp/src/comin.mod ]; then
    echo "whamp not built yet, building it now..."
    ( cd ../whamp/src && make ) || {
        echo "ERROR: whamp build failed (cd ../whamp/src && make)"; exit 1; }
fi

rm -f *.o *.mod                      # no stale modules from an earlier build
# modules first: everything else uses them
mpif90 $FLAGS $INC -c $(ls Module_Constants.f90 \
                           Module_Dispersion.f90 Module_Accounting.f90 2>/dev/null)
mpif90 $FLAGS *.f90 $INC -o "main_$name" $LIBS
echo "main_$name <- $(git rev-parse --short HEAD)$(git diff --quiet || echo ' +local changes')  flags: -O2 $*"
