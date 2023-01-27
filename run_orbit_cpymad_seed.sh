#!/bin/sh


mkdir tune_match_mb_fielderr_roll_mq_offset_90_IP5_10seeeds_corrhplus

for de in `seq 1 10`; do
    echo "$de" 
    sed -e "s/eseed=[0-9][0-9]/eseed=$de/g" test_cpymadx_tds.py > tmp.py
    sed -e "s/seed=[0-9][0-9];/seed=$de;/g" select_errors_tmp.madx > select_errors.madx
    python3 tmp.py > seed-"$de".out
    mv *.tab test_seed"$de"
    mv test_seed"$de"/ tune_match_mb_fielderr_roll_mq_offset_90_IP5_10seeeds_corrhplus/
    mv seed-"$de".out tune_match_mb_fielderr_roll_mq_offset_90_IP5_10seeeds_corrhplus/
    mv *.tab test_seed"$de"/
    rm tmp*.py select_errors.madx
done
