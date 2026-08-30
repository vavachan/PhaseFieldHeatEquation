#!/usr/bin/env bash

# file where dynamics information is dumped via SuperLattice::describe()
TARGET_IN=$1

# list containing C++ dynamics type information for performing cse
TARGET_OUT=$2

# define type of target
TYPE=dynamics

# Only look at optimizable dynamics with manageable complexity
awk -F\; '$3 == 1 && $6 > 0 && $6 < 100e6 { print $1 }' $TARGET_IN > tmp.$TYPE

# Perform some string manipulation of the following steps
awk 'NF > 0 { line=$0; gsub(/double/, "Expr", line);
              gsub(/float/, "Expr", line);
              line = gensub(/descriptors::D([23])Q([0-9]+)<[^>]+>/, "descriptors::D\\1Q\\2<>", "g", line);
              print line }' tmp.$TYPE | sort | uniq > "$TARGET_OUT"

rm tmp.$TYPE
