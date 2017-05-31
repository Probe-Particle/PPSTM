#!/bin/bash

echo "test for the PP-STM code:"
python PPSTM/PPAFM/generateLJFF.py -i TOAT.xyz --npy
python PPSTM/PPAFM/relaxed_scan.py --pos --npy
python PPSTM/PPAFM/plot_results.py --pos --df --save_df --npy
python PPSTM/dIdV_test_TOAT.py
echo "Now all things made, before submiting please run clean.sh"


