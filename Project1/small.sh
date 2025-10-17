#/usr/bin/bash

# for (( i=1; i<=100; i++ )); do
    julia --project=. src/run.jl data/small.csv 100 1000
# done