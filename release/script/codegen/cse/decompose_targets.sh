#!/usr/bin/env bash

# target file containing list of dynamics type information
EXT_TARGET_LIST=$1

# output directory
EXT_DIR=$2

# output for include header
INCLUDE_HEADER=$3

awk -v file="$INCLUDE_HEADER" -v dir="$EXT_DIR" 'BEGIN { print "/*  ========================================================" > file;
             print " *  ==  WARNING: This is an automatically generated file, ==" > file;
             print " *  ==                  do not modify.                    ==" > file;
             print " *  ========================================================" > file;
             print " */" > file;
             print "" > file;
             print "#ifndef DISABLE_CSE" > file}'

awk -v file="$INCLUDE_HEADER" -v dir="$EXT_DIR" '{ cmd = "echo -n \"" $0 "\" | sha256sum";
            cmd | getline hash;
            close(cmd);
            gsub(/ |-/, "", hash);
            print $0 > dir "/" hash ".txt";
            print "" >> file; print "//" $0 >> file; print "#include \"" hash ".cse.h\"" >> file
            }
            END { print "" > file; print "#endif" > file }' "$EXT_TARGET_LIST"

