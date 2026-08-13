"""
STM (Scanning Tunneling Microscopy) calculations using PyTorch.

This module provides PyTorch-based implementations for computing dI/dV (differential
conductance) spectra with sp(d)-sp(d) orbital interactions.
"""
from abc import ABC
from typing import Sequence, TypeAlias, List, Tuple

import numpy as np
import torch

from pyPPSTM.probe_stm_vectorized import ProbeStmVectorized, ProbeStmVectorizedSp, ProbeStmVectorizedSpd

DeviceLikeType: TypeAlias = str | torch.device | int


class ProbeStmPytorch(ProbeStmVectorized, ABC):
    """STM (Scanning Tunneling Microscopy) calculations using PyTorch."""
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
        """Compute dI/dV (differential conductance) using PyTorch.

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
            probe_stm = ProbeStmPytorchSpd(V, WF, eta, eig, R, Rat, coes, tip_coes)
        else:  # 4 sample orbitals = sp
            probe_stm = ProbeStmPytorchSp(V, WF, eta, eig, R, Rat, coes, tip_coes)
        return probe_stm().cpu()

    def __init__(self,
                 V: float,
                 WF: float,
                 eta: float,
                 eig: np.ndarray,
                 R: np.ndarray,
                 Rat: np.ndarray,
                 coes: np.ndarray,
                 tip_coes: np.ndarray):
        self._device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        super().__init__(V=V,
                        WF=WF,
                        eta=eta,
                        eig=torch.from_numpy(eig).to(self._device),
                        R=torch.from_numpy(R).to(self._device),
                        Rat=torch.from_numpy(Rat).to(self._device),
                        coes=torch.from_numpy(coes).to(self._device),
                        tip_coes=torch.from_numpy(tip_coes).to(self._device))

    @property
    def backend(self) -> str:
        return f"PyTorch ({self._device})"

    def _to_float(self, tensor):
        return tensor.float()

    def _unsqueeze(self, tensor, dims: int|List[int]|Tuple[int]):
        if isinstance(dims, int):
            tensor = torch.unsqueeze(tensor, dim=dims)
        else:
            new_shape = list(tensor.shape)
            for d in dims:
                new_shape.insert(d, 1)
            tensor = tensor.view(new_shape)
        return tensor

    def _zeros(self, size: Sequence[int], dtype):
        return torch.zeros(size, dtype=dtype, device=self._device)

    def _norm(self, tensor, ord: int, axis: int):
        return torch.linalg.norm(tensor, ord=ord, axis=axis)

    def _exp(self, tensor):
        return torch.exp(tensor)

    def _einsum(self, subscripts, *operands):
        return torch.einsum(subscripts, *operands)

    @property
    def _float32(self):
        return torch.float32


class ProbeStmPytorchSp(ProbeStmPytorch, ProbeStmVectorizedSp):
    pass


class ProbeStmPytorchSpd(ProbeStmPytorch, ProbeStmVectorizedSpd):
    pass
