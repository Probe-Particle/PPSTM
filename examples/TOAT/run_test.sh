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
python PPdos_simple.py
python ../../ppstm_run.py fixed.toml
python ../../ppstm_run.py relaxed_s.toml
python ../../ppstm_run.py relaxed_pxy.toml
echo "Now all things made, before submiting please run clean.sh"
