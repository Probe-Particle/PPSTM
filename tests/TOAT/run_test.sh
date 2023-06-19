#!/bin/bash
OMP=0	# 0 = 'False' , 1 = 'True'
if [ $OMP -eq 1 ]
then
    export OMP_NUM_THREADS=8
fi

echo "OMP_NUM_THREADS:"
echo $OMP_NUM_THREADS

echo "test for the PP-STM code:"
ppafm-generate-ljff -i TOAT.xyz -f npy
ppafm-relaxed-scan --pos -f npy
ppafm-plot-results --pos --df --save_df -f npy
python3 PPSTM/dIdV_test_TOAT.py
python3 PPdos_simple.py
python3 PPSTM_simple.py
echo "Now all things made, before submiting please run clean.sh"


