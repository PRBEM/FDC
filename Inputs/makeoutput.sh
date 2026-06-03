#!/bin/bash

SRC_DIR="/home/aroy/MySpace/Works/full-diffusion-code-f90-master/Inputs/LB_MLTs"
DEST_DIR="/home/aroy/MySpace/Works/full-diffusion-code-f90-master/Inputs/2026Chorus/UB"

mkdir -p "$DEST_DIR"

for i in $(seq 0 23); do
    mlttwo=$(printf "%02d" "$i")
    SRC_FILE="$SRC_DIR/params_MLT${i}.txt"
    DEST_FILE="$DEST_DIR/params_MLT${i}.txt"  # Use padded for destination

    if [[ ! -f "$SRC_FILE" ]]; then
        echo "Warning: $SRC_FILE not found, skipping."
        continue
    fi

    awk -v mlttwo="$mlttwo" '
    BEGIN { skip=0 }
    {
        # Identify the headers we want to replace and skip old values
        if ($0 ~ /^!!Edited by/) { print "!!Edited by Alwin Roy, Feb 26, 2026. If you have any question please contact Dedong Wang (dedongpotsdam@gmail.com)"; skip=0; next }
        if ($0 ~ /^!!!!!!!!!Output file name/) { print "!!!!!!!!!Output file name"; printf "AR_UB_MLT%s_Sheely01.plt\n", mlttwo; skip=1; next }
        if ($0 ~ /^!!!!!!!!!NPA,NE,NL/) { print "!!!!!!!!!NPA,NE,NL"; print "89 , 189 , 5"; skip=1; next }
        if ($0 ~ /^!!!!!!!!!L_min\*10,L_max\*10/) { print "!!!!!!!!!L_min*10,L_max*10"; print "30, 70"; skip=1; next }
        if ($0 ~ /^!!!!!!!!!E_min,Emax/) { print "!!!!!!!!!E_min,Emax in unit of keV"; print "0.001 , 30000"; skip=1; next }
        if ($0 ~ /^!!!!!!!!!DenModel,DenValue/) { print "!!!!!!!!!DenModel,DenValue in unit of cm^-3,DenLatINdex"; print "Sheely01 , 0.0 , 0.0"; skip=1; next }
        if ($0 ~ /^!!!!!!!!!Frequency Distribution/) { print "!!!!!!!!!Frequency Distribution"; print "2026UBChorus"; skip=1; next }
        if ($0 ~ /^!!!!!!!!!FreqUnit,omega_m_per,d_omega_per,o_lc,o_uc/) { print "!!!!!!!!!FreqUnit,omega_m_per,d_omega_per,o_lc,o_uc"; print "2026UBChorus, 0.2,  0.15,  0.05,  0.4"; skip=1; next }
        if ($0 ~ /^!!!!!!!!!lam_min,lam_max in unit of Degree/) { print "!!!!!!!!!lam_min,lam_max in unit of Degree"; print "0.0,  10.0"; skip=1; next }
        if ($0 ~ /^!!!!!!!!!BwType,Bw_Input/) { print "!!!!!!!!!BwType,Bw_Input in unit of pT"; print "2026UBChorus,  0.0"; skip=1; next }

        # Skip the old value line if it comes immediately after a header
        if (skip==1 && $0 ~ /^[0-9A-Za-z,. _-]+$/) { next }

        # Otherwise print the line
        print
        skip=0
    }
    ' "$SRC_FILE" > "$DEST_FILE"

    echo "Created modified file: $DEST_FILE"
done
