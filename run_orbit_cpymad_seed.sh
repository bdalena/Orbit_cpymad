#!/bin/sh

for de in `seq 1 10`; do
    echo "$de" 
    sed -e "s/eseed=[0-9][0-9]/eseed=$de/g" test_cpymadx_tds.py > tmp.py
    python3 tmp.py > seed-"$de".out
    mv *.tab test_seed"$de"/
    rm tmp*.py
done
