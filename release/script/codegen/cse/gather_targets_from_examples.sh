#!/usr/bin/env bash

find ../../../examples/ -name '*.dynamics' -exec sh -c '
for file; do
     parents=$(dirname "$file" | awk -F/ "{print \$(NF-2) \"_\" \$(NF-1) \"_\" \$(NF)}");
     cp "$file" "./targets/dynamics/${parents}_$(basename "$file")";
done
' sh {} +

find ../../../examples/ -name '*.operator' -exec sh -c '
for file; do
     parents=$(dirname "$file" | awk -F/ "{print \$(NF-2) \"_\" \$(NF-1) \"_\" \$(NF)}");
     cp "$file" "./targets/operator/${parents}_$(basename "$file")";
done
' sh {} +
