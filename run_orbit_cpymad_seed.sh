#!/bin/sh

mkdir mq_offset_20_IP5_2it_100seeeds_corrhplus

for de in `seq 1 100`; do
    echo "$de" 
    sed -e "s/eseed=[0-9][0-9]/eseed=$de/g" test_cpymadx_tds.py > tmp.py
    sed -e "s/seed=[0-9][0-9];/seed=$de;/g" select_errors_tmp.madx > select_errors.madx
    python3 tmp.py > seed-"$de".out
    mv *.tab test_seed"$de"
    mv test_seed"$de"/ mq_offset_20_IP5_2it_100seeeds_corrhplus/
    mv seed-"$de".out mq_offset_20_IP5_2it_100seeeds_corrhplus/
    rm tmp*.py select_errors.madx
done
