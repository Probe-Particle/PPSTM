#!/bin/bash

echo "test for the PP-STM code on the example of spin-polarized CuPc molecule, precalculated by FHI-AIMS:"
python PPSTM/dIdV_test_CuPc.py
python PPSTM/dIdV_test_CuPc_cp2k.py
echo "Now all things made, before submiting pleas run clean.sh"

