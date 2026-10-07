#!/bin/bash
# =============================================================================
# bisect_runs.sh - which commit changed the results?
#
# Builds every commit of the speed-up series in its own worktree, runs the same
# small input with each, and compares consecutive results. One short run per
# commit, so the whole chain takes minutes rather than days.
#
# Usage:  ./bisect_runs.sh /abs/path/to/Inputs/small_cold.txt [np]
#
# The input file must be an ABSOLUTE path and should be small
# (e.g. 23 pitch angles x 8 energies x 1 L shell). Its output file name is
# irrelevant: every run gets its own directory.
# =============================================================================
set -e

INPUT=${1:?usage: ./bisect_runs.sh /abs/path/to/input.txt [np]}
NP=${2:-28}
REPO=$(git rev-parse --show-toplevel)
# sources present in Src/ but not committed; copied into every worktree
UNTRACKED_SRC="Dump_WHAMP_Surface.f90 Module_Accounting.f90"
WORK=$HOME/bisect
CMP=$REPO/Src/compare_plt.py

# commit (or WORKING for the current, uncommitted state) + label
STEPS=(
  "17c8050:00_base_before_fixes"
  "f766044:01_cunningham_hot_fixes"
  "abd37c1:02_local_mode"
  "08e9a5b:03_profiling_timers"
  "807f172:04_root_finder_added"      # still off by default in this commit
  "c86709b:05_root_finder_default"
  "653e67a:06_system_clock"
  "136f0ad:07_vacuum_permittivity"
  "e37139d:08_fast_int_N"
  "63c8d9a:09_theta_skip"
  "WORKING:10_working_tree"           # G2 guard + mirror fix, uncommitted
)

mkdir -p "$WORK"

for step in "${STEPS[@]}"; do
    sha=${step%%:*}; label=${step##*:}
    src="$WORK/$label/Src"

    if [ ! -x "$WORK/$label/run/main_x" ]; then
        echo "=== building $label ($sha)"
        rm -rf "$WORK/$label"
        git worktree prune
        if [ "$sha" = "WORKING" ]; then
            mkdir -p "$WORK/$label"
            cp -r "$REPO/Src" "$WORK/$label/Src"
        else
            git worktree add --detach "$WORK/$label" "$sha" >/dev/null
        fi
        mkdir -p "$WORK/$label/run"
        # sources that are untracked in the repository would otherwise be missing
        # in the worktrees and the link step would fail
        for f in $UNTRACKED_SRC; do
            [ -f "$src/$f" ] || cp "$REPO/Src/$f" "$src/" 2>/dev/null || true
        done
        ( cd "$src" \
            && mpif90 -cpp -O2 -c $(ls Module_Constants.f90 Module_Variables.f90 \
                   Module_Dispersion.f90 Module_Accounting.f90 2>/dev/null) \
                   -I"$REPO"/whamp/src -ffree-line-length-none \
            && mpif90 -cpp -O2 *.f90 -I"$REPO"/whamp/src -o "$WORK/$label/run/main_x" \
                   -I/usr/lib/x86_64-linux-gnu/openmpi/include -ffree-line-length-none \
                   -L/usr/lib/ -l:libblas.so.3 -l:liblapack.so.3 "$REPO"/whamp/src/libwhamp.a \
              > "$WORK/$label/build.log" 2>&1 ) \
          || { echo "  BUILD FAILED - see $WORK/$label/build.log"; continue; }
    fi

    if [ -z "$(ls "$WORK/$label"/run/*.plt 2>/dev/null)" ]; then
        echo "=== running $label"
        ( cd "$WORK/$label/run" && /usr/bin/time -f "    %e s" \
              mpirun -np "$NP" ./main_x "$INPUT" > run.log 2>&1 ) || true
    fi
done

echo
echo "================ differences between consecutive steps ================"
prev=""
for step in "${STEPS[@]}"; do
    label=${step##*:}
    plt=$(ls "$WORK/$label"/run/*.plt 2>/dev/null | head -1)
    [ -z "$plt" ] && { echo "--- $label: no output"; continue; }
    if [ -n "$prev" ]; then
        echo
        echo "--- $prevlabel  ->  $label"
        python3 "$CMP" "$plt" "$prev" | tail -n +3
    fi
    prev=$plt; prevlabel=$label
done
