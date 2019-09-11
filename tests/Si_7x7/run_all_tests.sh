#!/bin/bash
OMP=	# 0 = 'False' , 1 = 'True'
if [ $OMP -eq 1 ]
then
    export OMP_NUM_THREADS=8
fi

echo "OMP_NUM_THREADS:"
echo $OMP_NUM_THREADS
echo "Now the tests:"

echo "test for the PP-STM code ploting DOS of Si(111) 7x7 reconstruction:"
python2 PPSTM/PDOS_test_Si_7x7.py
echo "test for the PP-STM code on the example os Si(111) 7x7 reconstruction - single STM calc -- +- 0.5 V:"
python2 PPSTM/dIdV_test_Si_7x7.py
echo "Now all things made, before submiting please run clean.sh"
echo "test for the PP-STM code on the example os Si(111) 7x7 reconstruction - full STM & dIdV Scan:"
python2 PPSTM/STM_test_Si_7x7.py
echo "Now all things made, before submiting please run clean.sh"


