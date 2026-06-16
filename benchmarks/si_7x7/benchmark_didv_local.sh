#!/bin/bash

echo "Run PP-STM on Si(111) 7x7 reconstruction - full STM & dIdV Scan:"

export OMP_NUM_THREADS=8
python ../../benchmark_didv.py fixed.toml exec_times_didv_si_7x7_fixed_cpp_8threads_local.csv
python ../../benchmark_didv.py fixed_wfd1.toml exec_times_didv_si_7x7_fixed_wfd1_cpp_8threads_local.csv

export OMP_NUM_THREADS=4
python ../../benchmark_didv.py fixed.toml exec_times_didv_si_7x7_fixed_cpp_4threads_local.csv
python ../../benchmark_didv.py fixed_wfd1.toml exec_times_didv_si_7x7_fixed_wfd1_cpp_4threads_local.csv

unset OMP_NUM_THREADS

python ../../benchmark_didv.py fixed_opencl.toml exec_times_didv_si_7x7_fixed_pyopencl_local.csv
python ../../benchmark_didv.py fixed_wfd1_opencl.toml exec_times_didv_si_7x7_fixed_wfd1_opencl_local.csv

python ../../benchmark_didv.py fixed_numpy.toml exec_times_didv_si_7x7_fixed_numpy_local.csv
python ../../benchmark_didv.py fixed_wfd1_numpy.toml exec_times_didv_si_7x7_fixed_wfd1_numpy_local.csv

python ../../benchmark_didv.py fixed_pytorch.toml exec_times_didv_si_7x7_fixed_pytorch_local.csv
python ../../benchmark_didv.py fixed_wfd1_pytorch.toml exec_times_didv_si_7x7_fixed_wfd1_pytorch_local.csv
