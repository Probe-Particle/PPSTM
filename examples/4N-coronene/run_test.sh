#!/bin/bash
OMP=0	# 0 = 'False' , 1 = 'True'
if [ $OMP -eq 1 ]
then
    export OMP_NUM_THREADS=8
fi

echo "OMP_NUM_THREADS:"
echo $OMP_NUM_THREADS
echo "Now the tests:"

which python3

echo "test for the PP-STM code:"
ppafm-generate-ljff -i crazy_mol.xyz -f npy
ppafm-relaxed-scan --pos -f npy
ppafm-plot-results --pos --df --save_df -f npy

python3 dIdV_test_4N-coronene.py
python3 dIdV_test_4N-coronene_cp2k.py
python3 ../../ppstm_run.py s-relaxed-png.toml
python3 ../../ppstm_run.py pxy-relaxed.toml
python3 ../../ppstm_run.py s-relaxed.toml
python3 PPdos_simple.py
python3 SUM_Orb_Contrib.py
echo "Now all things made, before submiting, please run clean.sh!"
