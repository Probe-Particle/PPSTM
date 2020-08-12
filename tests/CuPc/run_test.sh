#!/bin/bash
OMP=0	# 0 = 'False' , 1 = 'True'
if [ $OMP -eq 1 ]
then
    export OMP_NUM_THREADS=8
fi

echo "OMP_NUM_THREADS:"
echo $OMP_NUM_THREADS
echo "Now the tests:"

echo "test for the PP-STM code on the example of spin-polarized CuPc molecule, precalculated by FHI-AIMS & CP2K:"
python3 PPSTM/dIdV_test_CuPc.py
python3 PPSTM/dIdV_test_CuPc_cp2k.py
python3 PPdos_simple.py
python3 PPSTM_simple.py
echo "Now all things made, before submiting pleas run clean.sh"

