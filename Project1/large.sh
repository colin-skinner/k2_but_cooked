#/usr/bin/bash

# for (( i=1; i<=100; i++ )); do
  julia --project=. src/run.jl data/large.csv 1 300
# done