"""
STM (Scanning Tunneling Microscopy) calculations using Python.

This module provides Python-based implementations for computing dI/dV (differential
conductance) spectra with sp(d)-sp(d) orbital interactions.
"""
import logging
from abc import abstractmethod, ABC
from functools import partial
from typing import List, Callable, Sequence

import numpy as np
import math


class ProbeStmPython(ABC):
    """STM (Scanning Tunneling Microscopy) calculations using Python."""
    _SAMPLE_ORBITAL_COUNT: int
    _AB = 1.889725989
    _EV = 0.036749034
    _N_P = math.sqrt(3)
    _N_D = math.sqrt(15)
    _N_D2 = math.sqrt(5) * 0.5
    _I_3 = 0.3333333

    def __init__(self,
                 V: float,
                 WF: float,
                 eta: float,
                 eig,
                 R,
                 Rat,
                 coes,
                 tip_coes):
        self._v = np.float32(V)
        self._wf = np.float32(WF)
        self._eta = np.float32(eta)
        self._eig = self._to_float(eig)
        self._sample_atom_position = self._to_float(Rat)                       # shape                (n_a, 3)
        self._tip_position = self._unsqueeze(self._to_float(R), dim=3)         # shape (n_z, n_y, n_x,   1, 3)
        self._tip_coes = self._to_float(tip_coes)                              # shape (9,)
        self._coes = self._to_float(coes)\
                         .reshape(len(eig),
                                  len(Rat),
                                  self._SAMPLE_ORBITAL_COUNT)  # shape (n_o, n_a, orb_t)
        self._logger = logging.getLogger(self.__class__.__name__)

    def __call__(self):
        self._logger.debug(f"Entering the dI/dV ( sp(d)-sp(d) ) {self.backend} procedure")

        decay = math.sqrt((2 * self._wf) * self._EV)
        r_a = (self._sample_atom_position - self._tip_position) * self._AB             # shape (n_z, n_y, n_x, n_a, 3)
        r_a_norm = self._norm(r_a, ord=2, axis=-1)                                     # shape (n_z, n_y, n_x, n_a)

        radial_sq = self._exp(-decay * r_a_norm) ** 2                                  # shape (n_z, n_y, n_x, n_a)

        sample_dos = (self._eta /
            ((2 * math.pi) * ((self._eig - self._v) ** 2 + 0.25 * self._eta ** 2)))  # shape (n_o,), Lorentzian

        coes_sq_sample_dos = self._einsum("i...,i",
                                          self._coes ** 2,
                                          sample_dos)                                                 # shape (n_a, orb_t)

        g = self._zeros(self._tip_position.shape[:3], dtype=self._float32)
        for i, fn in enumerate(self.orbital_conductances_fns):
            if self._tip_coes[i] > 0.:
                orbital_conductances_i = fn(r_a, r_a_norm, coes_sq_sample_dos, decay)  # shape (n_z, n_y, n_x, n_a)

                g_atom = self._einsum("...a,...a",
                                      orbital_conductances_i,
                                      radial_sq)                                       # shape (n_z, n_y, n_x)

                g += self._tip_coes[i] * g_atom

        g *= 16 * math.pi ** 3 * decay

        self._logger.debug(f"dI/dV ( sp(d)-sp(d) ) {self.backend} procedure DONE")
        return g

    @property
    @abstractmethod
    def orbital_conductances_fns(self) -> List[Callable]:
        """Return a list of conductance functions for specific sample orbitals,
        where each function represents a different tip orbital.

        Returns:
            List[Callable]: functions computing orbital conductances.
                Each function must return arrays of shape (..., orb_t).
        """
        pass

    def _s(self, r_a: np.ndarray, r_a_norm: np.ndarray, coes_sq_sample_dos: np.ndarray, *args, **kwargs):
        """Orbital conductances for different sp or spd orbitals of sample on an s orbital tip.

        Args:
            r_a (torch.Tensor): Difference between tip position and sample atoms positions with shape (n_z, n_y, n_x, ..., 3).
            r_a_norm (torch.Tensor): Euclidean norm of r_a with shape (n_z, n_y, n_x, ...).
            coes_sq_sample_dos (torch.Tensor): shape (..., orb_t)

        Returns:
            torch.Tensor, shape (n_z, n_y, n_x, ...)
        """
        n_p_div_r_a_norm = self._N_P / r_a_norm
        n_d_div_r_a_norm_sq = self._N_D / r_a_norm ** 2

        t = coes_sq_sample_dos[:, 0]                                            # s  orb. of sample
        t = t + coes_sq_sample_dos[:, 1] * (n_p_div_r_a_norm * r_a[..., 1])**2  # py orb. of sample
        t = t + coes_sq_sample_dos[:, 2] * (n_p_div_r_a_norm * r_a[..., 2])**2  # pz orb. of sample
        t = t + coes_sq_sample_dos[:, 3] * (n_p_div_r_a_norm * r_a[..., 0])**2  # px orb. of sample
        if self._SAMPLE_ORBITAL_COUNT == 9:
            t = t + coes_sq_sample_dos[:, 4] * (n_d_div_r_a_norm_sq * r_a[..., 0] * r_a[..., 1])**2                   # dxy orb. of sample
            t = t + coes_sq_sample_dos[:, 5] * (n_d_div_r_a_norm_sq * r_a[..., 1] * r_a[..., 2])**2                   # dyz orb. of sample
            t = t + coes_sq_sample_dos[:, 6] * (self._N_D2 * (3 * r_a[..., 2] ** 2 / r_a_norm ** 2 - 1))**2           # dz2 orb. of sample
            t = t + coes_sq_sample_dos[:, 7] * (n_d_div_r_a_norm_sq * r_a[..., 0] * r_a[..., 2])**2                   # dxz orb. of sample
            t = t + coes_sq_sample_dos[:, 8] * (n_d_div_r_a_norm_sq * .5 * (r_a[..., 0] ** 2 - r_a[..., 1] ** 2))**2  # dx2-y2 orb. of sample
        return t

    def _py(self, r_a: np.ndarray, r_a_norm: np.ndarray, coes_sq_sample_dos: np.ndarray, decay: float):
        """Orbital conductances for different sp or spd orbitals of sample on a py orbital tip.

        Args:
            r_a (torch.Tensor): Difference between tip position and sample atoms positions with shape (n_z, n_y, n_x, ..., 3).
            r_a_norm (torch.Tensor): Euclidean norm of r_a with shape (n_z, n_y, n_x, ...).
            coes_sq_sample_dos (torch.Tensor): shape (..., orb_t)
            decay (float)

        Returns:
            torch.Tensor, shape (n_z, n_y, n_x, ...)
        """
        r_a_norm_inv = 1 / r_a_norm
        r_a_norm_inv_sq = r_a_norm_inv**2
        r_a_y_sq = r_a[..., 1]**2
        r_a_z_sq = r_a[..., 2]**2
        t =  coes_sq_sample_dos[:, 0] * (r_a[..., 1] * decay)**2                                                             # s  orb. of sample
        t += coes_sq_sample_dos[:, 1] * (self._N_P * ( -1 + decay*r_a_norm_inv*r_a_y_sq + r_a_norm_inv_sq*r_a_y_sq ))**2     # py orb. of sample
        t += coes_sq_sample_dos[:, 2] * (self._N_P * r_a[..., 2] * r_a[..., 2] * (decay*r_a_norm_inv + r_a_norm_inv_sq))**2  # pz orb. of sample
        t += coes_sq_sample_dos[:, 3] * (self._N_P * r_a[..., 2] * r_a[..., 0]*( decay*r_a_norm_inv + r_a_norm_inv_sq ))**2  # px orb. of sample
        if self._SAMPLE_ORBITAL_COUNT == 9:
            t += coes_sq_sample_dos[:, 4] * (self._N_D * r_a[..., 1] * r_a_norm_inv*( 2*r_a_y_sq*r_a_norm_inv_sq + decay*r_a_y_sq*r_a_norm_inv - 1 ))**2                # dxy orb. of sample
            t += coes_sq_sample_dos[:, 5] * (self._N_D * r_a[..., 2] * r_a_norm_inv*( 2*r_a_y_sq*r_a_norm_inv_sq + decay*r_a_y_sq*r_a_norm_inv - 1 ))**2                # dyz orb. of sample
            t += coes_sq_sample_dos[:, 6] * (self._N_D2 * ( (6*r_a[..., 1]**3*r_a_norm_inv**3-6*r_a[..., 2]*r_a_norm_inv) + decay*(3*r_a_z_sq*r_a_norm_inv_sq-1) ))**2  # dz2 orb. of sample
            t += coes_sq_sample_dos[:, 7] * (self._N_D * r_a[..., 0] * r_a[..., 1]*r_a[..., 2]*r_a_norm_inv_sq*( 2*r_a_norm_inv + decay ))**2                           # dxz orb. of sample
            t += coes_sq_sample_dos[:, 8] * (self._N_D * r_a_norm_inv_sq*r_a[..., 1]*((r_a[..., 0]**2-r_a_y_sq)*r_a_norm_inv_sq + 0.5*decay*(r_a[..., 0]**2-r_a_y_sq)*r_a_norm_inv + 1))**2  # dx2-y2 orb. of sample
        return t * r_a_norm_inv_sq

    def _pz(self, r_a: np.ndarray, r_a_norm: np.ndarray, coes_sq_sample_dos: np.ndarray, decay: float):
        """Orbital conductances for different sp or spd orbitals of sample on a pz orbital tip.

        Args:
            r_a (torch.Tensor): Difference between tip position and sample atoms positions with shape (n_z, n_y, n_x, ..., 3).
            r_a_norm (torch.Tensor): Euclidean norm of r_a with shape (n_z, n_y, n_x, ...).
            coes_sq_sample_dos (torch.Tensor): shape (..., orb_t)
            decay (float)

        Returns:
            torch.Tensor, shape (n_z, n_y, n_x, ...)
        """
        r_a_norm_inv = 1 / r_a_norm
        r_a_norm_inv_sq = r_a_norm_inv**2
        r_a_z_sq = r_a[..., 2]**2
        t =  coes_sq_sample_dos[:, 0] * (r_a[..., 2] * decay)**2                                                             # s  orb. of sample
        t += coes_sq_sample_dos[:, 1] * (self._N_P * r_a[..., 2] * r_a[..., 1] * (decay*r_a_norm_inv + r_a_norm_inv_sq))**2  # py orb. of sample
        t += coes_sq_sample_dos[:, 2] * (self._N_P * ( -1 + decay*r_a_norm_inv*r_a_z_sq + r_a_norm_inv_sq*r_a_z_sq ))**2     # pz orb. of sample
        t += coes_sq_sample_dos[:, 3] * (self._N_P * r_a[..., 2] * r_a[..., 0]*( decay*r_a_norm_inv + r_a_norm_inv_sq ))**2  # px orb. of sample
        if self._SAMPLE_ORBITAL_COUNT == 9:
            t += coes_sq_sample_dos[:, 4] * (self._N_D * r_a[..., 0] * r_a[..., 1]*r_a[..., 2]*r_a_norm_inv_sq*( 2*r_a_norm_inv + decay ))**2                           # dxy orb. of sample
            t += coes_sq_sample_dos[:, 5] * (self._N_D * r_a[..., 1] * r_a_norm_inv*( 2*r_a_z_sq*r_a_norm_inv_sq + decay*r_a_z_sq*r_a_norm_inv - 1 ))**2                # dyz orb. of sample
            t += coes_sq_sample_dos[:, 6] * (self._N_D2 * ( (6*r_a[..., 2]**3*r_a_norm_inv**3-6*r_a[..., 2]*r_a_norm_inv) + decay*(3*r_a_z_sq*r_a_norm_inv_sq-1) ))**2  # dz2 orb. of sample
            t += coes_sq_sample_dos[:, 7] * (self._N_D * r_a[..., 0] * r_a_norm_inv*( 2*r_a_z_sq*r_a_norm_inv_sq + decay*r_a_z_sq*r_a_norm_inv - 1 ))**2                # dxz orb. of sample
            t += coes_sq_sample_dos[:, 8] * (self._N_D * r_a_norm_inv_sq*r_a[..., 2]*(r_a[..., 0]**2-r_a[..., 1]**2)*( r_a_norm_inv + 0.5*decay  ))**2                  # dx2-y2 orb. of sample
        return t * r_a_norm_inv_sq

    def _px(self, r_a: np.ndarray, r_a_norm: np.ndarray, coes_sq_sample_dos: np.ndarray, decay: float):
        """Orbital conductances for different sp or spd orbitals of sample on a px orbital tip.

        Args:
            r_a (torch.Tensor): Difference between tip position and sample atoms positions with shape (n_z, n_y, n_x, ..., 3).
            r_a_norm (torch.Tensor): Euclidean norm of r_a with shape (n_z, n_y, n_x, ...).
            coes_sq_sample_dos (torch.Tensor): shape (..., orb_t)
            decay (float)

        Returns:
            torch.Tensor, shape (n_z, n_y, n_x, ...)
        """
        r_a_norm_inv = 1 / r_a_norm
        r_a_norm_inv_sq = r_a_norm_inv**2
        r_a_x_sq = r_a[..., 0]**2
        r_a_z_sq = r_a[..., 2]**2
        t =  coes_sq_sample_dos[:, 0] * (r_a[..., 0] * decay)**2                                                             # s  orb. of sample
        t += coes_sq_sample_dos[:, 1] * (self._N_P * r_a[..., 2] * r_a[..., 1] * (decay*r_a_norm_inv + r_a_norm_inv_sq))**2  # py orb. of sample
        t += coes_sq_sample_dos[:, 2] * (self._N_P * r_a[..., 2] * r_a[..., 2] * (decay*r_a_norm_inv + r_a_norm_inv_sq))**2  # pz orb. of sample
        t += coes_sq_sample_dos[:, 3] * (self._N_P * ( -1 + decay*r_a_norm_inv*r_a_x_sq + r_a_norm_inv_sq*r_a_x_sq ))**2     # px orb. of sample
        if self._SAMPLE_ORBITAL_COUNT == 9:
            t += coes_sq_sample_dos[:, 4] * (self._N_D * r_a[..., 1] * r_a_norm_inv*( 2*r_a_x_sq*r_a_norm_inv_sq + decay*r_a_x_sq*r_a_norm_inv - 1 ))**2                                     # dxy orb. of sample
            t += coes_sq_sample_dos[:, 5] * (self._N_D * r_a[..., 0] * r_a[..., 1]*r_a[..., 2]*r_a_norm_inv_sq*( 2*r_a_norm_inv + decay ))**2                                                # dyz orb. of sample
            t += coes_sq_sample_dos[:, 6] * (self._N_D2 * ( (6*r_a[..., 0]**3*r_a_norm_inv**3-6*r_a[..., 2]*r_a_norm_inv) + decay*(3*r_a_z_sq*r_a_norm_inv_sq-1) ))**2                       # dz2 orb. of sample
            t += coes_sq_sample_dos[:, 7] * (self._N_D * r_a[..., 2] * r_a_norm_inv*( 2*r_a_x_sq*r_a_norm_inv_sq + decay*r_a_x_sq*r_a_norm_inv - 1 ))**2                                     # dxz orb. of sample
            t += coes_sq_sample_dos[:, 8] * (self._N_D * r_a_norm_inv_sq*r_a[..., 0]*((r_a_x_sq-r_a[..., 1]**2)*r_a_norm_inv_sq + 0.5*decay*(r_a_x_sq-r_a[..., 1]**2)*r_a_norm_inv + 1))**2  # dx2-y2 orb. of sample
        return t * r_a_norm_inv_sq

    def _missing_conductance(self, tip_orbital: str, *args, **kwargs):
        """Placeholder returning zeros for unsupported tip-sample orbital combinations.

        Returns:
            torch.Tensor: Array of zeros, shape (n_z, n_y, n_x, n_a)
        """
        self._logger.info(f"Not implemented formula for {tip_orbital} tip orb"
                          f" and {'spd' if self._SAMPLE_ORBITAL_COUNT == 9 else 'sp'} sample orbitals!")
        return self._zeros(self._tip_position.shape[:-1] + len(self._sample_atom_position), dtype=self._float32)

    @property
    @abstractmethod
    def backend(self) -> str:
        pass

    @abstractmethod
    def _to_float(self, tensor):
        pass

    @abstractmethod
    def _unsqueeze(self, tensor, dim: int):
        pass

    @abstractmethod
    def _zeros(self, size: Sequence[int], dtype):
        pass

    @abstractmethod
    def _norm(self, tensor, ord: int, axis: int):
        pass

    @abstractmethod
    def _exp(self, tensor):
        pass

    @abstractmethod
    def _einsum(self, subscripts, *operands):
        pass

    @property
    @abstractmethod
    def _float32(self):
        pass


class ProbeStmPythonSp(ProbeStmPython, ABC):
    _SAMPLE_ORBITAL_COUNT = 4

    @property
    def orbital_conductances_fns(self) -> List[Callable]:
        """Return a list of conductance functions for sp sample orbitals,
        where each function represents a different tip orbital.

        Returns:
            List[Callable]: functions computing orbital conductances.
                Each function must return arrays of shape (..., orb_t).
        """
        return [
            self._s, self._py, self._pz, self._px,
            partial(self._missing_conductance, tip_orbital="dxy"),
            self._dyz, self._dz2, self._dxz,
            partial(self._missing_conductance, tip_orbital="dx2-y2"),
        ]

    def _dxz(self, r_a: np.ndarray, r_a_norm: np.ndarray, coes_sq_sample_dos: np.ndarray, decay: float):
        """Orbital conductances for different sp orbitals of sample on a dxz orbital tip.

        Args:
            r_a (torch.Tensor): Difference between tip position and sample atoms positions with shape (n_z, n_y, n_x, ..., 3).
            r_a_norm (torch.Tensor): Euclidean norm of r_a with shape (n_z, n_y, n_x, ...).
            coes_sq_sample_dos (torch.Tensor): shape (..., orb_t)
            decay (float)

        Returns:
            torch.Tensor, shape (n_z, n_y, n_x, ...)
        """
        decay_sq = decay**2
        r_a_norm_inv = 1 / r_a_norm
        r_a_norm_inv_sq = r_a_norm_inv**2
        r_a_norm_inv_cu = r_a_norm_inv**3
        r_a_x_sq = r_a[..., 0]**2
        r_a_z_sq = r_a[..., 2]**2
        r_a_x_mul_z = r_a[..., 0]*r_a[..., 2]
        t =  coes_sq_sample_dos[:, 0] * (r_a_x_mul_z*(r_a_norm_inv + decay))**2*decay_sq  # s  orb. of sample
        t += coes_sq_sample_dos[:, 1] * (self._N_P*r_a[..., 1]*r_a_x_mul_z*r_a_norm_inv*( 3*r_a_norm_inv_sq + 3*decay*r_a_norm_inv + decay_sq ))**2                                           # px orb. of sample
        t += coes_sq_sample_dos[:, 2] * (self._N_P*r_a[..., 0]*( 3*r_a_z_sq*r_a_norm_inv_cu + 3*decay*r_a_z_sq*r_a_norm_inv_sq - r_a_norm_inv + decay_sq*r_a_z_sq*r_a_norm_inv - decay ))**2  # pz orb. of sample
        t += coes_sq_sample_dos[:, 3] * (self._N_P*r_a[..., 2]*( 3*r_a_x_sq*r_a_norm_inv_cu + 3*decay*r_a_x_sq*r_a_norm_inv_sq - r_a_norm_inv + decay_sq*r_a_x_sq*r_a_norm_inv - decay ))**2  # py orb. of sample
        t *= r_a_norm_inv_sq**2
        return t

    def _dyz(self, r_a: np.ndarray, r_a_norm: np.ndarray, coes_sq_sample_dos: np.ndarray, decay: float):
        """Orbital conductances for different sp orbitals of sample on a dyz orbital tip.

        Args:
            r_a (torch.Tensor): Difference between tip position and sample atoms positions with shape (n_z, n_y, n_x, ..., 3).
            r_a_norm (torch.Tensor): Euclidean norm of r_a with shape (n_z, n_y, n_x, ...).
            coes_sq_sample_dos (torch.Tensor): shape (..., orb_t)
            decay (float)

        Returns:
            torch.Tensor, shape (n_z, n_y, n_x, ...)
        """
        decay_sq = decay**2
        r_a_norm_inv = 1 / r_a_norm
        r_a_norm_inv_sq = r_a_norm_inv**2
        r_a_norm_inv_cu = r_a_norm_inv**3
        r_a_y_sq = r_a[..., 1]**2
        r_a_z_sq = r_a[..., 2]**2
        r_a_y_mul_z = r_a[..., 1]*r_a[..., 2]
        t =  coes_sq_sample_dos[:, 0] * (r_a_y_mul_z*(r_a_norm_inv + decay))**2*decay_sq  # s  orb. of sample
        t += coes_sq_sample_dos[:, 1] * (self._N_P*r_a[..., 2]*( 3*r_a_y_sq*r_a_norm_inv_cu + 3*decay*r_a_y_sq*r_a_norm_inv_sq - r_a_norm_inv + decay_sq*r_a_y_sq*r_a_norm_inv - decay ))**2  # py orb. of sample
        t += coes_sq_sample_dos[:, 2] * (self._N_P*r_a[..., 1]*( 3*r_a_z_sq*r_a_norm_inv_cu + 3*decay*r_a_z_sq*r_a_norm_inv_sq - r_a_norm_inv + decay_sq*r_a_z_sq*r_a_norm_inv - decay ))**2  # pz orb. of sample
        t += coes_sq_sample_dos[:, 3] * (self._N_P*r_a[..., 0]*r_a_y_mul_z*r_a_norm_inv*( 3*r_a_norm_inv_sq + 3*decay*r_a_norm_inv + decay_sq ))**2  # px orb. of sample
        return t * r_a_norm_inv_sq**2

    def _dz2(self, r_a: np.ndarray, r_a_norm: np.ndarray, coes_sq_sample_dos: np.ndarray, decay: float):
        """Orbital conductances for different sp orbitals of sample on a dz2 orbital tip.

        Args:
            r_a (torch.Tensor): Difference between tip position and sample atoms positions with shape (n_z, n_y, n_x, ..., 3).
            r_a_norm (torch.Tensor): Euclidean norm of r_a with shape (n_z, n_y, n_x, ...).
            coes_sq_sample_dos (torch.Tensor): shape (..., orb_t)
            decay (float)

        Returns:
            torch.Tensor, shape (n_z, n_y, n_x, ...)
        """
        n_p_div_r_a_norm = self._N_P / r_a_norm
        decay_div_r_a_norm = decay / r_a_norm
        decay_sq = decay**2
        r_a_norm_inv = 1 / r_a_norm
        r_a_norm_inv_sq = r_a_norm_inv**2
        r_a_norm_inv_cu = r_a_norm_inv**3
        r_a_norm_inv_qu = r_a_norm_inv**4
        r_a_z_sq = r_a[..., 2]**2
        t =  coes_sq_sample_dos[:, 0] * (-self._I_3 * decay_sq + decay * r_a_z_sq * r_a_norm_inv_cu + decay_sq * r_a_z_sq * r_a_norm_inv_sq - decay_div_r_a_norm)**2  # s  orb. of sample
        t += coes_sq_sample_dos[:, 1] * (n_p_div_r_a_norm*r_a[..., 1]*( 3*r_a_z_sq*r_a_norm_inv_qu + 3*decay*r_a_z_sq*r_a_norm_inv_cu - r_a_norm_inv_sq + decay_sq*r_a_z_sq*r_a_norm_inv_sq - decay_div_r_a_norm - self._I_3*decay_sq))**2      # py orb. of sample
        t += coes_sq_sample_dos[:, 1] * (n_p_div_r_a_norm*r_a[..., 2]*( 3*r_a_z_sq*r_a_norm_inv_qu + 3*decay*r_a_z_sq*r_a_norm_inv_cu - 3*r_a_norm_inv_sq + decay_sq*r_a_z_sq*r_a_norm_inv_sq - 3*decay_div_r_a_norm - self._I_3*decay_sq))**2  # pz orb. of sample
        t += coes_sq_sample_dos[:, 1] * (n_p_div_r_a_norm*r_a[..., 0]*( 3*r_a_z_sq*r_a_norm_inv_qu + 3*decay*r_a_z_sq*r_a_norm_inv_cu - r_a_norm_inv_sq + decay_sq*r_a_z_sq*r_a_norm_inv_sq - decay_div_r_a_norm - self._I_3*decay_sq))**2      # px orb. of sample
        return t


class ProbeStmPythonSpd(ProbeStmPython, ABC):
    _SAMPLE_ORBITAL_COUNT = 9

    @property
    def orbital_conductances_fns(self) -> List[Callable]:
        """Return a list of conductance functions for spd sample orbitals,
        where each function represents a different tip orbital.

        Returns:
            List[Callable]: functions computing orbital conductances.
                Each function must return arrays of shape (..., orb_t).
        """
        return [
            self._s,
            self._py,
            self._pz,
            self._px,
            partial(self._missing_conductance, tip_orbital="dxy"),
            partial(self._missing_conductance, tip_orbital="dyz"),
            partial(self._missing_conductance, tip_orbital="dz2"),
            partial(self._missing_conductance, tip_orbital="dxz"),
            partial(self._missing_conductance, tip_orbital="dx2-y2"),
        ]
