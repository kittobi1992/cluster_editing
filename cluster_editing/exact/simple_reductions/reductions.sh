#!/bin/bash

# Compile the cpp converter
g++ -o reductions reductions.cpp

#for f in graphs/*.gr; do
#  ./reductions $f
#done

./reductions graphs/exact061.gr
