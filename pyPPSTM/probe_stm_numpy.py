"""
STM (Scanning Tunneling Microscopy) calculations using NumPy.

This module provides NumPy-based implementations for computing dI/dV (differential
conductance) spectra with sp(d)-sp(d) orbital interactions.
"""
from abc import ABC
from typing import Sequence, List, Tuple

import numpy as np

from pyPPSTM.probe_stm_python import ProbeStmPython, ProbeStmPythonSp, ProbeStmPythonSpd


class ProbeStmNumpy(ProbeStmPython, ABC):
    """STM (Scanning Tunneling Microscopy) calculations using NumPy."""
    _BACKEND = "NumPy"

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
             orb_t: int) -> np.ndarray:
        """Compute dI/dV (differential conductance) using NumPy.

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
        """
        if orb_t == 9:  # 9 sample orbitals = spd
            probe_stm = ProbeStmNumpySpd(V, WF, eta, eig, R, Rat, coes, tip_coes)
        else:  # 4 sample orbitals = sp
            probe_stm = ProbeStmNumpySp(V, WF, eta, eig, R, Rat, coes, tip_coes)
        return probe_stm()

    @property
    def backend(self) -> str:
        return self._BACKEND

    def _to_float(self, tensor):
        return tensor.astype(np.float32)

    def _unsqueeze(self, tensor, dims: int|List[int]|Tuple[int]):
        return np.expand_dims(tensor, axis=dims)

    def _zeros(self, size: Sequence[int], dtype):
        return np.zeros(size, dtype=dtype)

    def _norm(self, tensor, ord: int, axis: int):
        return np.linalg.norm(tensor, ord=ord, axis=axis)

    def _exp(self, tensor):
        return np.exp(tensor)

    def _einsum(self, subscripts, *operands):
        return np.einsum(subscripts, *operands, dtype=np.float32)

    @property
    def _float32(self):
        return np.float32


class ProbeStmNumpySp(ProbeStmNumpy, ProbeStmPythonSp):
    pass


class ProbeStmNumpySpd(ProbeStmNumpy, ProbeStmPythonSpd):
    pass
