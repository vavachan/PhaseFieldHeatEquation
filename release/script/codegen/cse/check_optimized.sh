#!/usr/bin/env bash

# location of the expression tree files
EXP_DIR=$1

# loaction of the optimized files
CSE_DIR=$2

COUNTER=0
COUNTER_FAILED=0
COUNTER_SUCCESS=0

for target_file in $EXP_DIR/*.out; do
    stripped_name="${target_file#$EXP_DIR/}" # Remove path prefix
    base_name="${stripped_name%.out}" # Remove ".txt" suffix
    out_file="$CSE_DIR/$base_name.cse.h" # Rebuild the path for obj_file
    let COUNTER++
    if [ ! -f "$out_file" ]; then
        echo "Failed to optimize: $target_file"
        let COUNTER_FAILED++
    fi
done

declare -i COUNTER_SUCCESS=COUNTER-COUNTER_FAILED 
echo "Summary: SUCCESS=$COUNTER_SUCCESS / FAILED=$COUNTER_FAILED / TOTAL=$COUNTER"
