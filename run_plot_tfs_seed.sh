#!/bin/sh

for de in `seq 1 100`; do
    echo "$de"
    sed -e "s/eseed=[0-9][0-9][0-9]/eseed=$de/g" plot_tfs.py > tmp.py
    python3 tmp.py
    rm tmp*.py
done
