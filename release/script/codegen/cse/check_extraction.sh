#!/usr/bin/env bash

# location of the target files
EXT_DIR=$1

# loaction of the generated expression trees
EXP_DIR=$2

COUNTER=0
COUNTER_FAILED=0
COUNTER_SUCCESS=0

for target_file in $EXT_DIR/*.txt; do
    stripped_name="${target_file#$EXT_DIR/}" # Remove path prefix
    base_name="${stripped_name%.txt}" # Remove ".txt" suffix
    out_file="$EXP_DIR/$base_name.out" # Rebuild the path for obj_file
    let COUNTER++
    if [ ! -f "$out_file" ]; then
        echo "Failed to generate: $target_file"
        let COUNTER_FAILED++
    fi
done

declare -i COUNTER_SUCCESS=COUNTER-COUNTER_FAILED 
echo "Summary: SUCCESS=$COUNTER_SUCCESS / FAILED=$COUNTER_FAILED / TOTAL=$COUNTER"
