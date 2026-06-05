#!/bin/bash
OMP=1	# 0 = 'False' , 1 = 'True'
if [ $OMP -eq 1 ]
then
    export OMP_NUM_THREADS=4
fi

echo "OMP_NUM_THREADS:"
echo $OMP_NUM_THREADS
echo "Now the tests:"

python3 ../../ppstm_run.py orbitals.toml
python3 ../../ppstm_run.py pxy_spd.toml
python3 ../../ppstm_run.py pz_spd.toml
python3 ../../ppstm_run.py dz2_sp.toml
python3 ../../ppstm_run.py dxyz_sp.toml


echo "Now all things made, before submiting, please run clean.sh!"
