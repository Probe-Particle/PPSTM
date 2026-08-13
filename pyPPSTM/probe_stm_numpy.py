"""
STM (Scanning Tunneling Microscopy) calculations using NumPy.

This module provides NumPy-based implementations for computing dI/dV (differential
conductance) spectra with sp(d)-sp(d) orbital interactions.
"""
from abc import ABC
from typing import Sequence

import numpy as np

from pyPPSTM.probe_stm_vectorized import ProbeStmVectorized, ProbeStmVectorizedSp, ProbeStmVectorizedSpd


class ProbeStmNumpy(ProbeStmVectorized, ABC):
    """STM (Scanning Tunneling Microscopy) calculations using NumPy."""
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
        return "NumPy"

    def _to_float(self, array: np.ndarray):
        return array.astype(np.float32)

    def _unsqueeze(self, array: np.ndarray, axis: int | Sequence[int]):
        return np.expand_dims(array, axis=axis)

    def _float_zeros(self, size: Sequence[int]):
        return np.zeros(size, dtype=np.float32)

    def _norm(self, array: np.ndarray, ord: int, axis: int):
        return np.linalg.norm(array, ord=ord, axis=axis)

    def _exp(self, array: np.ndarray):
        return np.exp(array)

    def _einsum(self, subscripts, *operands):
        return np.einsum(subscripts, *operands, dtype=np.float32)


class ProbeStmNumpySp(ProbeStmNumpy, ProbeStmVectorizedSp):
    pass


class ProbeStmNumpySpd(ProbeStmNumpy, ProbeStmVectorizedSpd):
    pass
