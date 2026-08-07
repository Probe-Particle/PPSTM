import logging
import sys
import time
from pathlib import Path

import pandas as pd
import tomli

import ppstm_run
from pyPPSTM import ProbeSTM, DidvBackend
from pyPPSTM.probe_stm_opencl import ProbeSTMOpenCLParallel
from pyPPSTM.probe_stm_numpy import ProbeStmNumpy
from pyPPSTM.probe_stm_pytorch import ProbeStmPytorch

if __name__ == "__main__":
     logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)

     # Get config file from command line
     config_file = sys.argv[1]

     output_csv_path = Path(sys.argv[2])

     # Load config file
     with open(config_file, 'rb') as f:
         config = tomli.load(f)

     print(f"Loaded config from {config_file}")
     print(f"Config: {config}")

     execution_times = []

     didv_backend = config['advanced']['didv_backend'] if 'didv_backend' in config['advanced'] else ProbeSTM.DidvBackend.CPP

     def didv_sp_sp_timed_wrapper(*args, **kwargs):
          start_time = time.perf_counter()
          result = didv_official(*args, **kwargs)
          end_time = time.perf_counter()
          diff_time = end_time - start_time
          logging.info(f"[BENCHMARK] {didv_official.__name__} took: {diff_time:.6f}s")
          execution_times.append(diff_time)
          return result

     if didv_backend == DidvBackend.OPENCL:
          ProbeSTMOpenCLParallel.get_instance()
          didv_official = ProbeSTMOpenCLParallel.didv
          ProbeSTMOpenCLParallel.didv = didv_sp_sp_timed_wrapper
     elif didv_backend == DidvBackend.CPP:
          didv_official = ProbeSTM.dIdV_sp_sp
          ProbeSTM.dIdV_sp_sp = didv_sp_sp_timed_wrapper
     elif didv_backend == DidvBackend.PYTORCH:
          didv_official = ProbeStmPytorch.didv
          ProbeStmPytorch.didv = didv_sp_sp_timed_wrapper
     else:  # didv_backend == DidvBackend.NUMPY
          didv_official = ProbeStmNumpy.didv
          ProbeStmNumpy.didv = didv_sp_sp_timed_wrapper

     ppstm_run.run_simulation(config)

     output_csv_path.parent.mkdir(parents=False, exist_ok=True)
     pd.Series(execution_times, name="exec_time_s").to_csv(output_csv_path, index=False, mode="w")
