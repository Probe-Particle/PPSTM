#!/bin/bash
OMP=1	# 0 = 'False' , 1 = 'True'
if [ $OMP -eq 1 ]
then
    export OMP_NUM_THREADS=4
fi

echo "OMP_NUM_THREADS:"
echo $OMP_NUM_THREADS
echo "Now the tests:"

ppstm-run orbitals.toml # s vs spd
ppstm-run s_sp.toml
ppstm-run pxy_sp.toml
ppstm-run pxy_spd.toml
ppstm-run pz_sp.toml
ppstm-run pz_spd.toml
ppstm-run dz2_sp.toml
ppstm-run dxyz_sp.toml
ppstm-run s_spy_high_eta.toml
ppstm-run high_wf.toml

echo "Now all things made, before submiting, please run clean.sh!"
