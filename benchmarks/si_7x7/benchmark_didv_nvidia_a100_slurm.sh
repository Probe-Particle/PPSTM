#!/bin/bash -l
#SBATCH --time=00:20:00
#SBATCH --gpus=a100:1
#SBATCH --partition=gpu-a100-80g
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1-1
#SBATCH --ntasks=1
#SBATCH --mem=1G
#SBATCH --output=benchmark_didv_nvidia_a100.out

module purge
source activate ppstm-dev

echo "Run PP-STM on Si(111) 7x7 reconstruction - full STM & dIdV Scan:"
python ../../benchmark_didv.py fixed_opencl.toml exec_times_didv_si_7x7_fixed_opencl_nvidia_a100.csv
python ../../benchmark_didv.py fixed_wfd1_opencl.toml exec_times_didv_si_7x7_fixed_wfd1_opencl_nvidia_a100.csv
python ../../benchmark_didv.py fixed_pytorch.toml exec_times_didv_si_7x7_fixed_pytorch_nvidia_a100.csv
python ../../benchmark_didv.py fixed_wfd1_pytorch.toml exec_times_didv_si_7x7_fixed_wfd1_pytorch_nvidia_a100.csv

