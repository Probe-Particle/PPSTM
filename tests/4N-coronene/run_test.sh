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
python3 PPSTM/PPAFM/generateLJFF.py -i crazy_mol.xyz --npy
python3 PPSTM/PPAFM/relaxed_scan.py --pos --npy
python3 PPSTM/PPAFM/plot_results.py --pos --df --save_df --npy
python3 PPSTM/dIdV_test_4N-coronene.py
python3 PPSTM/dIdV_test_4N-coronene_cp2k.py
python3 PPdos_simple.py
python3 PPSTM_simple.py
python3 PPSTM_simple_pxy.py
python3 PPSTM_simple_s.py
python3 SUM_Orb_Contrib.py
echo "Now all things made, before submiting, please run clean.sh!"
