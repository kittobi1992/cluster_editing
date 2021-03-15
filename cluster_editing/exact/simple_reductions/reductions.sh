#!/bin/bash

g++ -o reductions reductions.cpp -D_GLIBCXX_DEBUG

for f in graphs/*.gr; do
 ./reductions $f
done

# ./reductions graphs/exact041.gr
