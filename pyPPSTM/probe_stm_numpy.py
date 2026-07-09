"""
STM (Scanning Tunneling Microscopy) calculations using Numpy.

This module provides Numpy-based implementations for computing dI/dV (differential
conductance) spectra with sp(d)-sp(d) orbital interactions.
"""
import logging
from abc import abstractmethod, ABC
from functools import partial
from typing import List, Callable

import numpy as np


class ProbeStmNumpy(ABC):
    """STM (Scanning Tunneling Microscopy) calculations using Numpy."""
    _SAMPLE_ORBITAL_COUNT: int
    _AB = 1.889725989
    _EV = 0.036749034
    _N_P = np.sqrt(3)
    _N_D = np.sqrt(15)
    _N_D2 = np.sqrt(5) * .5
    _I_3 = 0.3333333

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
        """Compute dI/dV (differential conductance) using Numpy.

        Args:
            V (float): Applied voltage bias
            WF (float): Work function
            eta (float): Broadening parameter (energy smearing)
            eig (np.ndarray): Eigenvalues array, shape (n_o,)
            R (np.ndarray): Spatial grid of probe positions, shape (n_x, n_y, n_z, 3)
            Rat (np.ndarray): Atomic positions array, shape (n_a, 3)
            coes (np.ndarray): Orbital coefficients for sample, shape (n_o, n_a*orb_t)
            tip_coes (np.ndarray): Orbital coefficients for tip, shape (9,)
            orb_t (int): Orbital type identifier (4 or 9) for sample

        Returns:
            np.ndarray: dI/dV spectrum, shape (n_x, n_y, n_z)
        """
        if orb_t == 9:  # 9 sample orbitals = spd
            probe_stm = ProbeStmNumpySpd(V, WF, eta, eig, R, Rat, coes, tip_coes)
        else:  # 4 sample orbitals = sp
            probe_stm = ProbeStmNumpySp(V, WF, eta, eig, R, Rat, coes, tip_coes)
        return probe_stm()

    def __init__(self, V, WF, eta, eig: np.ndarray, R: np.ndarray, Rat: np.ndarray, coes: np.ndarray, tip_coes):
        self._v = np.float32(V)
        self._wf = np.float32(WF)
        self._eta = np.float32(eta)
        self._eig = eig.astype(np.float32)                                        # shape (n_o,)
        self._sample_atom_position = Rat.astype(np.float32)                       # shape                (n_a, 3)
        self._tip_position = np.expand_dims(R.astype(np.float32), 3)         # shape (n_x, n_y, n_z,   1, 3)
        self._tip_coes = tip_coes.astype(np.float32)                              # shape (9,)
        self._coes = coes.astype(np.float32)\
                         .reshape(len(eig),
                                  len(Rat),
                                  self._SAMPLE_ORBITAL_COUNT)  # shape (n_o, n_a, orb_t)
        self._logger = logging.getLogger(self.__class__.__name__)

    def __call__(self) -> np.ndarray:
        self._logger.debug("Entering the dI/dV ( sp(d)-sp(d) ) NumPy procedure")

        decay = np.sqrt((2 * self._wf) * self._EV)
        r_a = (self._sample_atom_position - self._tip_position) * self._AB       # shape (n_x, n_y, n_z, n_a, 3)
        r_a_norm = np.linalg.norm(r_a, ord=2, axis=-1)                           # shape (n_x, n_y, n_z, n_a)
        # r_a_norm = np.clip(r_a_norm, 1e-10, None)  # Prevent division by zero

        radial_sq = np.exp(-decay * r_a_norm) ** 2                               # shape (n_x, n_y, n_z, n_a)

        sample_dos = self._eta / ((2 * np.pi) * ((self._eig - self._v) ** 2 + 0.25 * self._eta ** 2))  # shape (n_o,), Lorentzian

        coes_sq_sample_dos = np.einsum("i...,i",
                                       self._coes ** 2,
                                       sample_dos)                                # shape (n_a, orb_t)

        g = np.zeros(self._tip_position.shape[:3], dtype=np.float32)
        for i, fn in enumerate(self.orbital_conductances_fns):
            if self._tip_coes[i] > 0.:
                orbital_conductances_i = fn(r_a, r_a_norm, decay)                          # shape (n_x, n_y, n_z, n_a, orb_t)
                orbital_conductances_i = np.einsum("...i,...i",
                                                   orbital_conductances_i ** 2,
                                                   coes_sq_sample_dos)                     # shape (n_x, n_y, n_z, n_a)

                g_atom = np.einsum("...a,...a",
                                   orbital_conductances_i,
                                   radial_sq)                                              # shape (n_x, n_y, n_z)

                g += self._tip_coes[i] * g_atom

        g *= 16 * np.pi ** 3 * decay

        self._logger.debug("dI/dV ( sp(d)-sp(d) ) NumPy procedure DONE")
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

    def _missing_conductance(self, tip_orbital: str, *args, **kwargs) -> np.ndarray:
        """Placeholder returning zeros for unsupported tip-sample orbital combinations.

        Returns:
            np.ndarray: Array of zeros, shape (..., orb_t)
        """
        self._logger.info(f"Not implemented formula for {tip_orbital} tip orb"
                          f" and {'spd' if self._SAMPLE_ORBITAL_COUNT == 9 else 'sp'} sample orbitals!")
        return np.zeros(self._tip_position.shape[:-1] + (self._SAMPLE_ORBITAL_COUNT,))

    def _s(self, r_a: np.ndarray, r_a_norm: np.ndarray, *args, **kwargs) -> np.ndarray:
        '''Orbital conductances for different sp or spd orbitals of sample on an s orbital tip.

        Args:
            r_a (np.ndarray): Difference between tip position and sample atoms positions with shape (..., 3).
            r_a_norm (np.ndarray): Euclidean norm of r_a with shape (n_x, n_y, n_z, ...).

        Returns:
            np.ndarray, shape (..., orb_t)
        '''
        n_p_div_r_a_norm = self._N_P / r_a_norm
        n_d_div_r_a_norm_sq = self._N_D / r_a_norm ** 2
        t = [
            np.ones(r_a.shape[:-1]),         # s  orb. of sample
            n_p_div_r_a_norm * r_a[..., 1],  # py orb. of sample
            n_p_div_r_a_norm * r_a[..., 2],  # pz orb. of sample
            n_p_div_r_a_norm * r_a[..., 0],  # px orb. of sample
        ]
        if self._SAMPLE_ORBITAL_COUNT == 9:
            t += [
                n_d_div_r_a_norm_sq * r_a[..., 0] * r_a[..., 1],                   # dxy orb. of sample
                n_d_div_r_a_norm_sq * r_a[..., 1] * r_a[..., 2],                   # dyz orb. of sample
                self._N_D2 * (3 * r_a[..., 2] ** 2 / r_a_norm ** 2 - 1),           # dz2 orb. of sample
                n_d_div_r_a_norm_sq * r_a[..., 0] * r_a[..., 2],                   # dxz orb. of sample
                n_d_div_r_a_norm_sq * .5 * (r_a[..., 0] ** 2 - r_a[..., 1] ** 2),  # dx2-y2 orb. of sample
            ]
        return np.stack(t, axis=-1, dtype=np.float32)


    def _py(self, r_a: np.ndarray, r_a_norm: np.ndarray, decay: float) -> np.ndarray:
        '''Orbital conductances for different sp or spd orbitals of sample on a py orbital tip.

        Args:
            r_a (np.ndarray): Difference between tip position and sample atoms positions with shape (..., 3).
            r_a_norm (np.ndarray): Euclidean norm of r_a with shape (n_x, n_y, n_z, ...).
            decay (float)

        Returns:
            np.ndarray, shape (..., orb_t)
        '''
        r_a_norm_inv = 1 / r_a_norm
        r_a_norm_inv_sq = r_a_norm_inv**2
        r_a_y_sq = r_a[..., 1]**2
        r_a_z_sq = r_a[..., 2]**2
        t = [
            r_a[..., 1] * decay,                                                             # s  orb. of sample
            self._N_P * ( -1 + decay*r_a_norm_inv*r_a_y_sq + r_a_norm_inv_sq*r_a_y_sq ),     # py orb. of sample
            self._N_P * r_a[..., 2] * r_a[..., 2] * (decay*r_a_norm_inv + r_a_norm_inv_sq),  # pz orb. of sample
            self._N_P * r_a[..., 2] * r_a[..., 0]*( decay*r_a_norm_inv + r_a_norm_inv_sq ),  # px orb. of sample
        ]
        if self._SAMPLE_ORBITAL_COUNT == 9:
            t += [
                self._N_D * r_a[..., 1] * r_a_norm_inv*( 2*r_a_y_sq*r_a_norm_inv_sq + decay*r_a_y_sq*r_a_norm_inv - 1 ),                # dxy orb. of sample
                self._N_D * r_a[..., 2] * r_a_norm_inv*( 2*r_a_y_sq*r_a_norm_inv_sq + decay*r_a_y_sq*r_a_norm_inv - 1 ),                # dyz orb. of sample
                self._N_D2 * ( (6*r_a[..., 1]**3*r_a_norm_inv**3-6*r_a[..., 2]*r_a_norm_inv) + decay*(3*r_a_z_sq*r_a_norm_inv_sq-1) ),  # dz2 orb. of sample
                self._N_D * r_a[..., 0] * r_a[..., 1]*r_a[..., 2]*r_a_norm_inv_sq*( 2*r_a_norm_inv + decay ),                           # dxz orb. of sample
                self._N_D * r_a_norm_inv_sq*r_a[..., 1]*((r_a[..., 0]**2-r_a_y_sq)*r_a_norm_inv_sq + 0.5*decay*(r_a[..., 0]**2-r_a_y_sq)*r_a_norm_inv + 1),  # dx2-y2 orb. of sample
            ]
        return np.stack(t, axis=-1, dtype=np.float32) * np.expand_dims(r_a_norm_inv, axis=-1)


    def _pz(self, r_a: np.ndarray, r_a_norm: np.ndarray, decay: float) -> np.ndarray:
        '''Orbital conductances for different sp or spd orbitals of sample on a pz orbital tip.

        Args:
            r_a (np.ndarray): Difference between tip position and sample atoms positions with shape (..., 3).
            r_a_norm (np.ndarray): Euclidean norm of r_a with shape (n_x, n_y, n_z, ...).
            decay (float)

        Returns:
            np.ndarray, shape (..., orb_t)
        '''
        r_a_norm_inv = 1 / r_a_norm
        r_a_norm_inv_sq = r_a_norm_inv**2
        r_a_z_sq = r_a[..., 2]**2
        t = [
            r_a[..., 2] * decay,                                                             # s  orb. of sample
            self._N_P * r_a[..., 2] * r_a[..., 1] * (decay*r_a_norm_inv + r_a_norm_inv_sq),  # py orb. of sample
            self._N_P * ( -1 + decay*r_a_norm_inv*r_a_z_sq + r_a_norm_inv_sq*r_a_z_sq ),     # pz orb. of sample
            self._N_P * r_a[..., 2] * r_a[..., 0]*( decay*r_a_norm_inv + r_a_norm_inv_sq ),  # px orb. of sample
        ]
        if self._SAMPLE_ORBITAL_COUNT == 9:
            t += [
                self._N_D * r_a[..., 0] * r_a[..., 1]*r_a[..., 2]*r_a_norm_inv_sq*( 2*r_a_norm_inv + decay ),                           # dxy orb. of sample
                self._N_D * r_a[..., 1] * r_a_norm_inv*( 2*r_a_z_sq*r_a_norm_inv_sq + decay*r_a_z_sq*r_a_norm_inv - 1 ),                # dyz orb. of sample
                self._N_D2 * ( (6*r_a[..., 2]**3*r_a_norm_inv**3-6*r_a[..., 2]*r_a_norm_inv) + decay*(3*r_a_z_sq*r_a_norm_inv_sq-1) ),  # dz2 orb. of sample
                self._N_D * r_a[..., 0] * r_a_norm_inv*( 2*r_a_z_sq*r_a_norm_inv_sq + decay*r_a_z_sq*r_a_norm_inv - 1 ),                # dxz orb. of sample
                self._N_D * r_a_norm_inv_sq*r_a[..., 2]*(r_a[..., 0]**2-r_a[..., 1]**2)*( r_a_norm_inv + 0.5*decay  ),                  # dx2-y2 orb. of sample
            ]
        return np.stack(t, axis=-1, dtype=np.float32) * np.expand_dims(r_a_norm_inv, axis=-1)


    def _px(self, r_a: np.ndarray, r_a_norm: np.ndarray, decay: float) -> np.ndarray:
        '''Orbital conductances for different sp or spd orbitals of sample on a px orbital tip.

        Args:
            r_a (np.ndarray): Difference between tip position and sample atoms positions with shape (..., 3).
            r_a_norm (np.ndarray): Euclidean norm of r_a with shape (n_x, n_y, n_z, ...).
            decay (float)

        Returns:
            np.ndarray, shape (..., orb_t)
        '''
        r_a_norm_inv = 1 / r_a_norm
        r_a_norm_inv_sq = r_a_norm_inv**2
        r_a_x_sq = r_a[..., 0]**2
        r_a_z_sq = r_a[..., 2]**2
        t = [
            r_a[..., 0] * decay,                                                             # s  orb. of sample
            self._N_P * r_a[..., 2] * r_a[..., 1] * (decay*r_a_norm_inv + r_a_norm_inv_sq),  # py orb. of sample
            self._N_P * r_a[..., 2] * r_a[..., 2] * (decay*r_a_norm_inv + r_a_norm_inv_sq),  # pz orb. of sample
            self._N_P * ( -1 + decay*r_a_norm_inv*r_a_x_sq + r_a_norm_inv_sq*r_a_x_sq ),     # px orb. of sample
        ]
        if self._SAMPLE_ORBITAL_COUNT == 9:
            t += [
                self._N_D * r_a[..., 1] * r_a_norm_inv*( 2*r_a_x_sq*r_a_norm_inv_sq + decay*r_a_x_sq*r_a_norm_inv - 1 ),                                     # dxy orb. of sample
                self._N_D * r_a[..., 0] * r_a[..., 1]*r_a[..., 2]*r_a_norm_inv_sq*( 2*r_a_norm_inv + decay ),                                                # dyz orb. of sample
                self._N_D2 * ( (6*r_a[..., 0]**3*r_a_norm_inv**3-6*r_a[..., 2]*r_a_norm_inv) + decay*(3*r_a_z_sq*r_a_norm_inv_sq-1) ),                       # dz2 orb. of sample
                self._N_D * r_a[..., 2] * r_a_norm_inv*( 2*r_a_x_sq*r_a_norm_inv_sq + decay*r_a_x_sq*r_a_norm_inv - 1 ),                                     # dxz orb. of sample
                self._N_D * r_a_norm_inv_sq*r_a[..., 0]*((r_a_x_sq-r_a[..., 1]**2)*r_a_norm_inv_sq + 0.5*decay*(r_a_x_sq-r_a[..., 1]**2)*r_a_norm_inv + 1),  # dx2-y2 orb. of sample
            ]
        return np.stack(t, axis=-1, dtype=np.float32) * np.expand_dims(r_a_norm_inv, axis=-1)


class ProbeStmNumpySp(ProbeStmNumpy):
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

    def _dxz(self, r_a: np.ndarray, r_a_norm: np.ndarray, decay: float) -> np.ndarray:
        '''Orbital conductances for different sp orbitals of sample on a dxz orbital tip.

        Args:
            r_a (np.ndarray): Difference between tip position and sample atoms positions with shape (..., 3).
            r_a_norm (np.ndarray): Euclidean norm of r_a with shape (n_x, n_y, n_z, ...).
            decay (float)

        Returns:
            np.ndarray, shape (..., orb_t)
        '''
        decay_sq = decay**2
        r_a_norm_inv = 1 / r_a_norm
        r_a_norm_inv_sq = r_a_norm_inv**2
        r_a_norm_inv_cu = r_a_norm_inv**3
        r_a_x_sq = r_a[..., 0]**2
        r_a_z_sq = r_a[..., 2]**2
        r_a_x_mul_z = r_a[..., 0]*r_a[..., 2]
        t = [
            r_a_x_mul_z*decay*(r_a_norm_inv + decay),  # s  orb. of sample
            self._N_P*r_a[..., 1]*r_a_x_mul_z*r_a_norm_inv*( 3*r_a_norm_inv_sq + 3*decay*r_a_norm_inv + decay_sq ),                                           # px orb. of sample
            self._N_P*r_a[..., 0]*( 3*r_a_z_sq*r_a_norm_inv_cu + 3*decay*r_a_z_sq*r_a_norm_inv_sq - r_a_norm_inv + decay_sq*r_a_z_sq*r_a_norm_inv - decay ),  # pz orb. of sample
            self._N_P*r_a[..., 2]*( 3*r_a_x_sq*r_a_norm_inv_cu + 3*decay*r_a_x_sq*r_a_norm_inv_sq - r_a_norm_inv + decay_sq*r_a_x_sq*r_a_norm_inv - decay ),  # py orb. of sample
        ]
        return np.stack(t, axis=-1, dtype=np.float32) * np.expand_dims(r_a_norm_inv_sq, axis=-1)


    def _dyz(self, r_a: np.ndarray, r_a_norm: np.ndarray, decay: float) -> np.ndarray:
        '''Orbital conductances for different sp orbitals of sample on a dyz orbital tip.

        Args:
            r_a (np.ndarray): Difference between tip position and sample atoms positions with shape (..., 3).
            r_a_norm (np.ndarray): Euclidean norm of r_a with shape (n_x, n_y, n_z, ...).
            decay (float)

        Returns:
            np.ndarray, shape (..., orb_t)
        '''
        decay_sq = decay**2
        r_a_norm_inv = 1 / r_a_norm
        r_a_norm_inv_sq = r_a_norm_inv**2
        r_a_norm_inv_cu = r_a_norm_inv**3
        r_a_y_sq = r_a[..., 1]**2
        r_a_z_sq = r_a[..., 2]**2
        r_a_y_mul_z = r_a[..., 1]*r_a[..., 2]
        t = [
            r_a_y_mul_z*decay*(r_a_norm_inv + decay),  # s  orb. of sample
            self._N_P*r_a[..., 2]*( 3*r_a_y_sq*r_a_norm_inv_cu + 3*decay*r_a_y_sq*r_a_norm_inv_sq - r_a_norm_inv + decay_sq*r_a_y_sq*r_a_norm_inv - decay ),  # py orb. of sample
            self._N_P*r_a[..., 1]*( 3*r_a_z_sq*r_a_norm_inv_cu + 3*decay*r_a_z_sq*r_a_norm_inv_sq - r_a_norm_inv + decay_sq*r_a_z_sq*r_a_norm_inv - decay ),  # pz orb. of sample
            self._N_P*r_a[..., 0]*r_a_y_mul_z*r_a_norm_inv*( 3*r_a_norm_inv_sq + 3*decay*r_a_norm_inv + decay_sq ),  # px orb. of sample
        ]
        return np.stack(t, axis=-1, dtype=np.float32) * np.expand_dims(r_a_norm_inv_sq, axis=-1)


    def _dz2(self, r_a: np.ndarray, r_a_norm: np.ndarray, decay: float) -> np.ndarray:
        '''Orbital conductances for different sp orbitals of sample on a dz2 orbital tip.

        Args:
            r_a (np.ndarray): Difference between tip position and sample atoms positions with shape (..., 3).
            r_a_norm (np.ndarray): Euclidean norm of r_a with shape (n_x, n_y, n_z, ...).
            decay (float)

        Returns:
            np.ndarray, shape (..., orb_t)
        '''
        n_p_div_r_a_norm = self._N_P / r_a_norm
        decay_div_r_a_norm = decay / r_a_norm
        decay_sq = decay**2
        r_a_norm_inv = 1 / r_a_norm
        r_a_norm_inv_sq = r_a_norm_inv**2
        r_a_norm_inv_cu = r_a_norm_inv**3
        r_a_norm_inv_qu = r_a_norm_inv**4
        r_a_z_sq = r_a[..., 2]**2
        t = [
            -self._I_3 * decay_sq + decay * r_a_z_sq * r_a_norm_inv_cu + decay_sq * r_a_z_sq * r_a_norm_inv_sq - decay_div_r_a_norm,  # s  orb. of sample
            n_p_div_r_a_norm*r_a[..., 1]*( 3*r_a_z_sq*r_a_norm_inv_qu + 3*decay*r_a_z_sq*r_a_norm_inv_cu - r_a_norm_inv_sq + decay_sq*r_a_z_sq*r_a_norm_inv_sq - decay_div_r_a_norm - self._I_3*decay_sq),      # py orb. of sample
            n_p_div_r_a_norm*r_a[..., 2]*( 3*r_a_z_sq*r_a_norm_inv_qu + 3*decay*r_a_z_sq*r_a_norm_inv_cu - 3*r_a_norm_inv_sq + decay_sq*r_a_z_sq*r_a_norm_inv_sq - 3*decay_div_r_a_norm - self._I_3*decay_sq),  # pz orb. of sample
            n_p_div_r_a_norm*r_a[..., 0]*( 3*r_a_z_sq*r_a_norm_inv_qu + 3*decay*r_a_z_sq*r_a_norm_inv_cu - r_a_norm_inv_sq + decay_sq*r_a_z_sq*r_a_norm_inv_sq - decay_div_r_a_norm - self._I_3*decay_sq),      # px orb. of sample
        ]
        return np.stack(t, axis=-1, dtype=np.float32)


class ProbeStmNumpySpd(ProbeStmNumpy):
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
