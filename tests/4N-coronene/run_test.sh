#!/bin/bash
OMP=0	# 0 = 'False' , 1 = 'True'
if [ $OMP -eq 1 ]
then
    export OMP_NUM_THREADS=8
fi

echo "OMP_NUM_THREADS:"
echo $OMP_NUM_THREADS
echo "Now the tests:"

echo "test for the PP-STM code:"
python PPSTM/PPAFM/generateLJFF.py -i crazy_mol.xyz --npy
python PPSTM/PPAFM/relaxed_scan.py --pos --npy
python PPSTM/PPAFM/plot_results.py --pos --df --save_df --npy
python PPSTM/dIdV_test_4N-coronene.py
python PPSTM/dIdV_test_4N-coronene_cp2k.py
echo "Now all things made, before submiting, please run clean.sh!"
