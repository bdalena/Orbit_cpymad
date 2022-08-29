#!/bin/sh

for de in `seq 1 100`; do
    echo "$de" 
    sed -e "s/eseed=[0-9][0-9]/eseed=$de/g" test_cpymadx_tds.py > tmp.py
    sed -e "s/seed=[0-9][0-9];/seed=$de;/g" select_errors_tmp.madx > select_errors.madx
    python3 tmp.py > seed-"$de".out
    mv *.tab test_seed"$de"/
    rm tmp*.py select_errors.madx
done
