#!/usr/bin/env bash

# Define the input path
TARGET_IN=$1

# Define the output file
TARGET_OUT=$2

# Define type of targets
TYPE=$3

# Clear the output file if it already exists
> "$TARGET_OUT"

# Loop through all .log files in the current directory
for file in "$TARGET_IN"/*.$TYPE; do
    if [ -f "$file" ]; then
        echo "$file"
        tail +2 "$file" >> "$TARGET_OUT"
    fi
done
