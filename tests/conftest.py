import pytest

from pyPPSTM.ProbeSTM import dIdV_sp_sp as dIdV_cpp
from pyPPSTM.probe_stm_opencl import ProbeSTMOpenCLParallel, ProbeSTMOpenCLSequential
from pyPPSTM.probe_stm_numpy import ProbeStmNumpy
from pyPPSTM.probe_stm_pytorch import ProbeStmPytorch


@pytest.fixture(params=[
    dIdV_cpp,
    ProbeSTMOpenCLParallel.get_instance().didv,
    ProbeSTMOpenCLSequential.get_instance().didv,
    ProbeStmNumpy.didv,
    ProbeStmPytorch.didv,
], ids=[
    "C++",
    "OpenCL parallel",
    "OpenCL sequential",
    "NumPy",
    "PyTorch"
])
def didv_backend(request):
    return request.param
