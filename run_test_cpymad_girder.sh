#!/bin/sh

args=("$@")
sd=${args[0]}
rm -rf /tmp/bdalena/wrk_"$sd"; mkdir -p /tmp/bdalena/wrk_"$sd"
cp /feynman/work/dacm/leda/bdalena/FCC-ee/orbit_corr/Orbit_cpymad/test_cpymadx_tds_girder.py /tmp/bdalena/wrk_"$sd"/.
cp /feynman/work/dacm/leda/bdalena/FCC-ee/orbit_corr/Orbit_cpymad/def_functions2.py /tmp/bdalena/wrk_"$sd"/.
cp /feynman/work/dacm/leda/bdalena/FCC-ee/orbit_corr/Orbit_cpymad/select_errors_girder.madx /tmp/bdalena/wrk_"$sd"/.
cd /tmp/bdalena/wrk_"$sd"
sed -e "s/eseed=[0-9][0-9]/eseed=$sd/g" test_cpymadx_tds_girder.py > tmp1.py
sed -e "s/seed=[0-9][0-9];/seed=$sd;/g" tmp1.py > tmp.py
python3 tmp.py > seed-"$sd".out
mv *.tab test_seed"$sd"/
mv seed-"$sd".out test_seed"$sd"/
mv test_seed"$sd"/ /feynman/work/dacm/leda/bdalena/FCC-ee/orbit_corr/Orbit_cpymad/

cd ..
rm -rf wrk_"$sd"
