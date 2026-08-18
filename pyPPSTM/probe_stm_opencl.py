"""
GPU-accelerated STM (Scanning Tunneling Microscopy) calculations using OpenCL.

This module provides OpenCL-based implementations for computing dI/dV (differential
conductance) spectra with sp(d)-sp(d) orbital interactions. It offers both parallel
and sequential execution modes for different hardware configurations.
"""
import logging
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Callable

import numpy as np
import pyopencl as cl
from pyopencl import Buffer


class _ProbeSTMOpenCL(ABC):
    """GPU acceleration for STM (Scanning Tunneling Microscopy) calculations using OpenCL.

    This is an abstract base class that manages OpenCL context, command queue, and kernel
    execution. Subclasses must implement the execution strategy via `didv_opencl_fn`
    and `_global_size()`.

    Attributes:
        _INSTANCE: Singleton instance (class variable)
        _LOCAL_SIZE: Work group size for OpenCL (class variable)
    """
    _INSTANCE = None
    _LOCAL_SIZE = (1,)

    @classmethod
    def didv(cls,
             V: float,
             WF: float,
             eta: float,
             eig: np.ndarray,
             R: np.ndarray,
             Rat: np.ndarray,
             coes: np.ndarray,
             tip_coes: np.ndarray,
             orb_t: int):
        """Compute dI/dV (differential conductance) using OpenCL sp(d)-sp(d) kernel.

        Args:
            V (float): Applied voltage bias
            WF (float): Work function
            eta (float): Broadening parameter (energy smearing)
            eig (np.ndarray): Eigenvalues array, shape (n_e,)
            R (np.ndarray): Spatial grid of probe positions, shape (n_z, n_y, n_x, 3)
            Rat (np.ndarray): Atomic positions array, shape (n_a, 3)
            coes (np.ndarray): Orbital coefficients for sample, shape (n_e, n_a*orb_t)
            tip_coes (np.ndarray): Orbital coefficients for tip, shape (9,)
            orb_t (int): Orbital type identifier (4 or 9) for sample

        Returns:
            np.ndarray: dI/dV spectrum, shape (n_z, n_y, n_x)

        Note:
            - Arrays are converted to C-contiguous format
            - See ProbeSTM_spd.cl for kernel implementation details
        """
        self = cls.get_instance()

        self._logger.debug("Entering the dI/dV ( sp(d)-sp(d) ) OpenCL procedure")

        n_points = int(np.prod(R.shape[:3]))
        didv_1d = np.zeros(n_points, dtype=np.float32)
        didv_1d_ocl = self._prepare_write_only_tensor_for_opencl(didv_1d)

        self.didv_opencl_fn(self._ocl_queue,
                            self._global_size(n_points),
                            self._LOCAL_SIZE,
                            *list(map(np.int32, [orb_t, len(Rat), len(eig), n_points])),
                            *list(map(np.float32, [V, WF, eta])),
                            *list(map(self._prepare_read_only_float_tensor_for_opencl, [eig, R, Rat, coes, tip_coes])),
                            didv_1d_ocl)

        cl.enqueue_copy(self._ocl_queue, didv_1d, didv_1d_ocl)

        didv_1d = didv_1d.reshape(*R.shape[:3])

        self._logger.debug("OpenCL kernel DONE, we're back in Python")
        return didv_1d

    @classmethod
    def get_instance(cls, platform: int = 0):
        if cls._INSTANCE is None:
            cls._INSTANCE = cls(platform)
        return cls._INSTANCE

    def __init__(self, platform: int = 0):
        self._ocl_environment = self._OCLEnvironment(platform)
        self._ocl_program = self._ocl_environment.loadProgram(self._ocl_environment.CL_PATH / "ProbeSTM_spd.cl")
        self._ocl_context = self._ocl_environment.context
        self._ocl_queue = self._ocl_environment.queue
        self._ocl_kernel = None
        self._logger = logging.getLogger(self.__class__.__name__)

    def _prepare_read_only_float_tensor_for_opencl(self, tensor) -> Buffer:
        """Create read-only buffer from host array.

        Ensures C-contiguous layout and float32 dtype, then creates
        a buffer with automatic data copy.

        Args:
            tensor (np.ndarray): Host array (any floating point precision/layout)

        Returns:
            cl.Buffer: buffer with READ_ONLY flag set
        """
        if tensor.flags["C_CONTIGUOUS"]:
            tensor = tensor.astype(np.float32)
        else:
            tensor = np.ascontiguousarray(tensor, dtype=np.float32)
        return cl.Buffer(hostbuf=tensor,
                         flags=cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR,
                         context=self._ocl_context)

    def _prepare_write_only_tensor_for_opencl(self, tensor) -> Buffer:
        """Create write-only buffer to host array.

        Args:
            tensor (np.ndarray): Host array (float32/1D)

        Returns:
            cl.Buffer: buffer with WRITE_ONLY flag set
        """
        return cl.Buffer(size=tensor.nbytes,
                         flags=cl.mem_flags.WRITE_ONLY,
                         context=self._ocl_context)

    @property
    def didv_opencl_fn(self) -> Callable:
        """Reference to the OpenCL kernel function for dI/dV computation."""
        return self._ocl_kernel

    @abstractmethod
    def _global_size(self, n_points: int):
        """Define OpenCL global work size for kernel launch.

        Args:
            n_points (int): Total number of points to process

        Returns:
            tuple: Global work size, shape (m,)
        """
        pass

    class _OCLEnvironment:
        """Manages OpenCL context, platform, and device initialization.

        Attributes:
            CL_PATH (Path): Directory containing .cl kernel source files
            platform: Selected OpenCL platform
            context: OpenCL context for memory and kernel management
            queue: Command queue for kernel execution
        """
        CL_PATH = Path(__file__).parent.parent.joinpath("cl")

        def __init__(self, i_platform=0):
            platforms = self._get_platforms()
            try:
                self._platform = platforms[i_platform]
            except IndexError:
                raise ValueError(f"Device {i_platform} not found. Max available platform id: {len(platforms)-1}")
            self._context = cl.Context(properties=[(cl.context_properties.PLATFORM, self._platform)], devices=None)
            self._queue = cl.CommandQueue(self._context)
            self._logger = logging.getLogger(self.__class__.__name__)

        def __del__(self):
            self._queue.finish()
            self._logger.debug("Finishing OCL queue...")

        @property
        def context(self):
            return self._context

        @property
        def queue(self):
            return self._queue

        def loadProgram(self, fname):
            cl_path = str(self.CL_PATH)
            if self._platform.name != "Portable Computing Language":
                cl_path = f'"{cl_path}"'
            with open(fname) as f:
                program = cl.Program(self._context, f.read()).build(options=["-I", cl_path])
            return program

        @staticmethod
        def _get_platforms():
            try:
                platforms = cl.get_platforms()
            except cl._cl.LogicError:
                raise RuntimeError(
                    "Could not find any OpenCL platforms. Check that the OpenCL ICD for your device is installed.")
            return platforms


class ProbeSTMOpenCLParallel(_ProbeSTMOpenCL):
    """Parallel execution strategy for dI/dV calculations.

    Launches work-items for up to n_points simultaneously.
    """
    def __init__(self, platform: int = 0):
        super().__init__(platform=platform)
        self._ocl_kernel = self._ocl_program.proc_didv_spd_spd

    def _global_size(self, n_points: int):
        return (n_points,)


class ProbeSTMOpenCLSequential(_ProbeSTMOpenCL):
    """Sequential execution strategy for dI/dV calculations.

    Launches single work-item that processes points serially.
    """
    def __init__(self, platform: int = 0):
        super().__init__(platform=platform)
        self._ocl_kernel = self._ocl_program.proc_didv_spd_spd_sequential

    def _global_size(self, n_points: int):
        return (1,)
