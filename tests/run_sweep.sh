#!/bin/bash
# ==============================================================================
# run_sweep.sh
#
# Runs the theta_max sweep (HOT and COLD, for each theta_max value) using the
# existing compiled main_chorus_test binary. Input files are expected to sit
# in ../Inputs/sweep/, named {hot,cold}_thetamax{VALUE}.txt, matching the
# "Output file name" field already set inside each input file
# (sweep_{hot,cold}_thetamax{VALUE}.plt), so no renaming is needed afterward.
#
# Runs sequentially -- each run already uses -np 28 processes, so running
# multiple sweep points in parallel would oversubscribe most machines. If you
# have significantly more than 28 cores and want concurrency instead, see the
# commented-out parallel variant near the bottom.
#
# Usage: bash run_sweep.sh
# (run from the Src/ directory, same place you'd normally invoke mpirun)
# ==============================================================================

set -e   # stop on first error, rather than silently continuing past a failed run

NPROC=28
INPUT_DIR="../Inputs/sweep"
LOG_DIR="../Outputs/sweep_logs"
THETA_MAX_VALUES="10 20 30 35 40 42 45"
WAVE_TYPES="hot cold"

mkdir -p "$LOG_DIR"

echo "Starting theta_max sweep: ${THETA_MAX_VALUES}"
echo "Log files will be written to ${LOG_DIR}/"
echo ""

start_time=$(date +%s)

for theta_max in $THETA_MAX_VALUES; do
    for type in $WAVE_TYPES; do
        input_file="${INPUT_DIR}/${type}_thetamax${theta_max}.txt"
        log_file="${LOG_DIR}/${type}_thetamax${theta_max}.log"

        if [ ! -f "$input_file" ]; then
            echo "WARNING: ${input_file} not found, skipping."
            continue
        fi

        echo "[$(date +%H:%M:%S)] Running ${type}_thetamax${theta_max} ..."
        mpirun -np $NPROC ./main_chorus_test "$input_file" > "$log_file" 2>&1

        if [ $? -eq 0 ]; then
            echo "  -> done."
        else
            echo "  -> FAILED (see ${log_file})"
        fi
    done
done

end_time=$(date +%s)
elapsed=$((end_time - start_time))
echo ""
echo "Sweep complete in ${elapsed}s. Logs in ${LOG_DIR}/"

# -----------------------------------------------------------------------------
# OPTIONAL: parallel variant, if you have enough cores (e.g. NPROC=28 x up to
# N simultaneous jobs). Uncomment and comment out the sequential loop above.
# Requires GNU parallel or just backgrounding with wait; example with wait:
#
# for theta_max in $THETA_MAX_VALUES; do
#     for type in $WAVE_TYPES; do
#         input_file="${INPUT_DIR}/${type}_thetamax${theta_max}.txt"
#         log_file="${LOG_DIR}/${type}_thetamax${theta_max}.log"
#         mpirun -np $NPROC ./main_chorus_test "$input_file" > "$log_file" 2>&1 &
#     done
# done
# wait
# -----------------------------------------------------------------------------
