#!/bin/bash
OMP=0	# 0 = 'False' , 1 = 'True'
if [ $OMP -eq 1 ]
then
    export OMP_NUM_THREADS=8
fi

echo "OMP_NUM_THREADS:"
echo $OMP_NUM_THREADS
echo "Now the tests:"

python PPdos_simple.py

echo "Run PP-STM on Si(111) 7x7 reconstruction - full STM & dIdV Scan:"
python ../../ppstm_run.py fixed.toml
python ../../ppstm_run.py fixed_wfd1.toml
echo "Now all things made, before submiting please run clean.sh"


