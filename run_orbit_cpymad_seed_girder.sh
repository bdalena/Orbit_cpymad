#!/bin/sh

for de in `seq 1 100`; do
    echo "seed = $de"
    sbatch -n 1 run_test_cpymad_girder.sh $de
done
