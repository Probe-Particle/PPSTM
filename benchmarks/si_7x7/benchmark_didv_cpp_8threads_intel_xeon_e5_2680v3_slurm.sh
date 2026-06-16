#!/bin/bash -l
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1-1
#SBATCH --ntasks=1
#SBATCH --mem=1536M
#SBATCH --partition=batch-hsw
#SBATCH --constraint="hsw avx2"
#SBATCH --output=benchmark_didv_cpp_8threads_intel_xeon_e5_2680v3.out

module purge
source activate ppstm-dev

export OMP_NUM_THREADS=8

echo "Run PP-STM on Si(111) 7x7 reconstruction - full STM & dIdV Scan:"
python ../../benchmark_didv.py fixed.toml exec_times_didv_si_7x7_fixed_cpp_8threads_intel_xeon_e5_2680v3.csv
python ../../benchmark_didv.py fixed_wfd1.toml exec_times_didv_si_7x7_fixed_wfd1_cpp_8threads_intel_xeon_e5_2680v3.csv
