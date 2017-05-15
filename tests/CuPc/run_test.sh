#!/bin/bash

echo "test for the PP-STM code on the example of spin-polarized CuPc molecule, precalculated by FHI-AIMS:"
python PPSTM/dIdV_test_CuPc.py
#gdb -ex r --args python PPSTM/dIdV_test_CuPc.py
echo "Now all things made"
echo "Done, bye bye"


