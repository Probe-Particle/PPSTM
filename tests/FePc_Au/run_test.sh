#!/bin/bash
OMP=0	# 0 = 'False' , 1 = 'True'
if [ $OMP -eq 1 ]
then
    export OMP_NUM_THREADS=8
fi

echo "OMP_NUM_THREADS:"
echo $OMP_NUM_THREADS
echo "Now the tests:"

echo "test for the PP-STM code: IETS"
python2 PPSTM/PPAFM/generateLJFF.py -i geom-cube.in --npy
python2 PPSTM/PPAFM/relaxed_scan.py --pos --npy --vib 3
python2 PPSTM/PPAFM/plot_results.py --pos --df --iets 16 0.0015 0.001 --npy --WSxM --save_df

python2 PPSTM/IETS_test_FePc.py

echo "Now all things made, before submiting, please run clean.sh!"
