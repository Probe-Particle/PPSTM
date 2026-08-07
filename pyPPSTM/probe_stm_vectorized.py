"""
STM (Scanning Tunneling Microscopy) calculations using vectorization.

This module provides vectorized implementations for computing dI/dV (differential
conductance) spectra with sp(d)-sp(d) orbital interactions.
"""
import logging
import math
from abc import abstractmethod, ABC
from functools import partial
from typing import List, Callable, Sequence, Tuple

import numpy as np


class ProbeStmVectorized(ABC):
    """STM (Scanning Tunneling Microscopy) calculations using vectorization."""
    _SAMPLE_ORBITAL_COUNT: int
    _AB = 1.889725989
    _EV = 0.036749034
    _N_P = math.sqrt(3)
    _N_D = math.sqrt(15)
    _N_D2 = math.sqrt(5) * 0.5
    _I_3 = 1/3

    def __init__(self,
                 V: float,
                 WF: float,
                 eta: float,
                 eig,
                 R,
                 Rat,
                 coes,
                 tip_coes,
                 n_tip_position_chunks: int = 1):
        """Args:
            V (float): Applied voltage bias
            WF (float): Work function
            eta (float): Broadening parameter (energy smearing)
            eig (array_like): Eigenvalues array, shape (n_e,)
            R (array_like): Spatial grid of probe positions, shape (n_z, n_y, n_x, 3)
            Rat (array_like): Atomic positions array, shape (n_a, 3)
            coes (array_like): Orbital coefficients for sample, shape (n_e, n_a*orb_t)
            tip_coes (array_like): Orbital coefficients for tip, shape (9)
            orb_t (int): Orbital type identifier (4 or 9) for sample
            n_tip_position_chunks (int): nr. subsets of tip positions, default 1
        """
        self._logger = logging.getLogger(self.__class__.__name__)

        v = np.float32(V)
        wf = np.float32(WF)
        eta = np.float32(eta)
        self._out_shape = R.shape[:3]
        n_zyx = np.prod(self._out_shape)
        self._tip_dos = self._to_float(tip_coes)  # shape (9,)
        eig = self._to_float(eig)                                                                # shape      (n_e,)
        self._sample_dos = self._lorentzian(x=eig, loc=v, scale=0.5 * eta)                       # shape      (n_e,)

        sample_atom_position = self._to_float(Rat)                                               # shape           (n_a, 3)
        tip_position = self._unsqueeze(self._to_float(R), axis=[3, 4]).reshape(n_zyx, 1, 1, 3)   # shape (n_c,   1,   1, 3)

        self._decay = math.sqrt((2 * wf) * self._EV)
        self._r_a = (tip_position - sample_atom_position) * self._AB                             # shape (n_c,   1, n_a, 3)
        self._r_a_norm = self._euclidean_norm(self._r_a, axis=-1)                                # shape (n_c,   1, n_a)
        self._radial = self._exp(-self._decay * self._r_a_norm)                                  # shape (n_c,   1, n_a)

        self._coes = self._to_float(coes)\
                         .reshape(len(eig),
                                  len(Rat),
                                  self._SAMPLE_ORBITAL_COUNT)                                    # shape      (n_e, n_a,    orb_t)

        self._n_tip_position_chunks = n_tip_position_chunks
        self._tip_position_chunk_size = math.ceil(self._r_a.shape[0] / n_tip_position_chunks)
        self._logger.debug(f"{self._n_tip_position_chunks} chunks => chunk size: {self._tip_position_chunk_size}")

    def __call__(self):
        """Compute dI/dV."""

        self._logger.debug(f"Entering the dI/dV ( sp(d)-sp(d) ) {self.backend} procedure")

        g_tip_orbs_sum = self._float_zeros((self._r_a.shape[0], self._coes.shape[0]))
        for i, g_tip_orb_fn in enumerate(self.g_tip_orb_fns):
            if self._tip_dos[i] > 0.:
                g_tip_orb = []
                for j, from_j in enumerate(range(0, self._r_a.shape[0], self._tip_position_chunk_size)):
                    to_j = from_j + self._tip_position_chunk_size
                    g_tip_orb.append(
                        self._broadcasted_dot_last_axis(
                            g_tip_orb_fn(from_j, to_j),                                  # shape (n_c, n_e, n_a)
                            self._radial[from_j:to_j]                                    # shape (n_c,   1, n_a)
                        )                                                                # shape (n_c, n_e)
                    )
                    self._logger.debug(f"Done {j+1}/{self._n_tip_position_chunks} chunks")

                g_tip_orbs_sum += self._tip_dos[i] * self._vstack(g_tip_orb) ** 2

        g = self._broadcasted_dot_last_axis(g_tip_orbs_sum, self._sample_dos)\
            * 16 * math.pi ** 3 * self._decay                                            # shape (n_c)

        self._logger.debug(f"dI/dV ( sp(d)-sp(d) ) {self.backend} procedure DONE")
        return g.reshape(self._out_shape)

    @property
    @abstractmethod
    def g_tip_orb_fns(self) -> List[Callable]:
        """Return a list of conductance functions for specific sample orbitals,
        where each function represents a different tip orbital.

        Returns:
            List[Callable]: functions computing orbital conductances.
                Each function must return arrays of shape (n_c, n_e, n_a).
        """
        pass

    def _s(self, from_i, to_i):
        """Orbital conductances for different sp or spd orbitals of sample on an s orbital tip.

        Returns:
            array_like, shape (n_c, n_e, n_a)
        """
        r_a = self._r_a[from_i:to_i]            # shape (n_c, 1, n_a, 3)
        r_a_norm = self._r_a_norm[from_i:to_i]  # shape (n_c, 1, n_a)
        n_p_div_r_a_norm = self._N_P / r_a_norm
        r_a_norm_sq = r_a_norm ** 2
        n_d_div_r_a_norm_sq = self._N_D / r_a_norm_sq

        t = self._coes[..., 0]                                       # s  orb. of sample, shape      (n_e, n_a)
        t = t + self._coes[..., 1] * n_p_div_r_a_norm * r_a[..., 1]  # py orb. of sample, shape (n_c, n_e, n_a)
        t = t + self._coes[..., 2] * n_p_div_r_a_norm * r_a[..., 2]  # pz orb. of sample
        t = t + self._coes[..., 3] * n_p_div_r_a_norm * r_a[..., 0]  # px orb. of sample
        if self._SAMPLE_ORBITAL_COUNT == 9:
            t = t + self._coes[..., 4] * n_d_div_r_a_norm_sq * r_a[..., 0] * r_a[..., 1]  # dxy orb. of sample
            t = t + self._coes[..., 5] * n_d_div_r_a_norm_sq * r_a[..., 1] * r_a[..., 2]  # dyz orb. of sample
            t = t + self._coes[..., 6] * self._N_D2 * (3 * r_a[..., 2] ** 2 / r_a_norm_sq - 1)  # dz2 orb. of sample
            t = t + self._coes[..., 7] * n_d_div_r_a_norm_sq * r_a[..., 0] * r_a[..., 2]  # dxz orb. of sample
            t = t + self._coes[..., 8] * n_d_div_r_a_norm_sq * .5 * (r_a[..., 0] ** 2 - r_a[..., 1] ** 2)  # dx2-y2 orb. of sample
        return t

    def _py(self, from_i, to_i):
        """Orbital conductances for different sp or spd orbitals of sample on a py orbital tip.

        Returns:
            array_like, shape (n_c, n_e, n_a)
        """
        r_a = self._r_a[from_i:to_i]  # shape (n_c, 1, n_a, 3)
        r_a_norm = self._r_a_norm[from_i:to_i]  # shape (n_c, 1, n_a)
        r_a_norm_inv = 1 / r_a_norm
        r_a_norm_inv_sq = r_a_norm_inv**2
        r_a_x_sq = r_a[..., 0]**2
        r_a_y_sq = r_a[..., 1]**2
        r_a_z_sq = r_a[..., 2]**2
        t =  self._coes[..., 0] * r_a[..., 1] * self._decay  # s  orb. of sample, shape (n_c, n_e, n_a)
        t += self._coes[..., 1] * self._N_P * ( -1 + self._decay*r_a_norm_inv*r_a_y_sq + r_a_norm_inv_sq*r_a_y_sq )     # py orb. of sample
        t += self._coes[..., 2] * self._N_P * r_a[..., 1] * r_a[..., 2] * (self._decay*r_a_norm_inv + r_a_norm_inv_sq)  # pz orb. of sample
        t += self._coes[..., 3] * self._N_P * r_a[..., 1] * r_a[..., 0] * (self._decay*r_a_norm_inv + r_a_norm_inv_sq)  # px orb. of sample
        if self._SAMPLE_ORBITAL_COUNT == 9:
            t += self._coes[..., 4] * self._N_D * r_a[..., 0] * r_a_norm_inv*( 2*r_a_y_sq*r_a_norm_inv_sq + self._decay*r_a_y_sq*r_a_norm_inv - 1 )   # dxy orb. of sample
            t += self._coes[..., 5] * self._N_D * r_a[..., 2] * r_a_norm_inv*( 2*r_a_y_sq*r_a_norm_inv_sq + self._decay*r_a_y_sq*r_a_norm_inv - 1 )   # dyz orb. of sample
            t += self._coes[..., 6] * self._N_D2 * r_a[..., 1] * (6*r_a_z_sq*r_a_norm_inv**3 + self._decay*(3*r_a_z_sq*r_a_norm_inv_sq-1))            # dz2 orb. of sample
            t += self._coes[..., 7] * self._N_D * r_a[..., 0] * r_a[..., 1]*r_a[..., 2]*r_a_norm_inv_sq*( 2*r_a_norm_inv + self._decay )  # dxz orb. of sample
            t += self._coes[..., 8] * self._N_D * r_a_norm_inv*r_a[..., 1]*((r_a_x_sq-r_a_y_sq)*r_a_norm_inv_sq + 0.5*self._decay*(r_a_x_sq-r_a_y_sq)*r_a_norm_inv + 1)  # dx2-y2 orb. of sample
        return t * r_a_norm_inv

    def _pz(self, from_i, to_i):
        """Orbital conductances for different sp or spd orbitals of sample on a pz orbital tip.

        Returns:
            array_like, shape (n_c, n_e, n_a)
        """
        r_a = self._r_a[from_i:to_i]  # shape (n_c, 1, n_a, 3)
        r_a_norm = self._r_a_norm[from_i:to_i]  # shape (n_c, 1, n_a)
        r_a_norm_inv = 1 / r_a_norm
        r_a_norm_inv_sq = r_a_norm_inv**2
        r_a_z_sq = r_a[..., 2]**2
        t =  self._coes[..., 0] * r_a[..., 2] * self._decay  # s  orb. of sample, shape (n_c, n_e, n_a)
        t += self._coes[..., 1] * self._N_P * r_a[..., 2] * r_a[..., 1] * (self._decay*r_a_norm_inv + r_a_norm_inv_sq)  # py orb. of sample
        t += self._coes[..., 2] * self._N_P * ( -1 + self._decay*r_a_norm_inv*r_a_z_sq + r_a_norm_inv_sq*r_a_z_sq )                 # pz orb. of sample
        t += self._coes[..., 3] * self._N_P * r_a[..., 2] * r_a[..., 0]*( self._decay*r_a_norm_inv + r_a_norm_inv_sq )  # px orb. of sample
        if self._SAMPLE_ORBITAL_COUNT == 9:
            t += self._coes[..., 4] * self._N_D * r_a[..., 0] * r_a[..., 1]*r_a[..., 2]*r_a_norm_inv_sq*( 2*r_a_norm_inv + self._decay )                   # dxy orb. of sample
            t += self._coes[..., 5] * self._N_D * r_a[..., 1] * r_a_norm_inv*( 2*r_a_z_sq*r_a_norm_inv_sq + self._decay*r_a_z_sq*r_a_norm_inv - 1 )                    # dyz orb. of sample
            t += self._coes[..., 6] * self._N_D2 * ((6*r_a[..., 2]**3*r_a_norm_inv**3-6*r_a[..., 2]*r_a_norm_inv) + self._decay*(3*r_a_z_sq*r_a_norm_inv_sq-1))  # dz2 orb. of sample
            t += self._coes[..., 7] * self._N_D * r_a[..., 0] * r_a_norm_inv*( 2*r_a_z_sq*r_a_norm_inv_sq + self._decay*r_a_z_sq*r_a_norm_inv - 1 )                    # dxz orb. of sample
            t += self._coes[..., 8] * self._N_D * r_a_norm_inv_sq*r_a[..., 2]*(r_a[..., 0]**2-r_a[..., 1]**2)*( r_a_norm_inv + 0.5*self._decay  )          # dx2-y2 orb. of sample
        return t * r_a_norm_inv

    def _px(self, from_i, to_i):
        """Orbital conductances for different sp or spd orbitals of sample on a px orbital tip.

        Returns:
            array_like, shape (n_c, n_e, n_a)
        """
        r_a = self._r_a[from_i:to_i]  # shape (n_c, 1, n_a, 3)
        r_a_norm = self._r_a_norm[from_i:to_i]  # shape (n_c, 1, n_a)
        r_a_norm_inv = 1 / r_a_norm
        r_a_norm_inv_sq = r_a_norm_inv**2
        r_a_x_sq = r_a[..., 0]**2
        r_a_y_sq = r_a[..., 1]**2
        r_a_z_sq = r_a[..., 2]**2
        t =  self._coes[..., 0] * r_a[..., 0] * self._decay  # s  orb. of sample, shape (n_c, n_e, n_a)
        t += self._coes[..., 1] * self._N_P * r_a[..., 0] * r_a[..., 1] * (self._decay*r_a_norm_inv + r_a_norm_inv_sq)  # py orb. of sample
        t += self._coes[..., 2] * self._N_P * r_a[..., 0] * r_a[..., 2] * (self._decay*r_a_norm_inv + r_a_norm_inv_sq)  # pz orb. of sample
        t += self._coes[..., 3] * self._N_P * ( -1 + self._decay*r_a_norm_inv*r_a_x_sq + r_a_norm_inv_sq*r_a_x_sq )                 # px orb. of sample
        if self._SAMPLE_ORBITAL_COUNT == 9:
            t += self._coes[..., 4] * self._N_D * r_a[..., 1] * r_a_norm_inv*( 2*r_a_x_sq*r_a_norm_inv_sq + self._decay*r_a_x_sq*r_a_norm_inv - 1 )                      # dxy orb. of sample
            t += self._coes[..., 5] * self._N_D * r_a[..., 0] * r_a[..., 1]*r_a[..., 2]*r_a_norm_inv_sq*( 2*r_a_norm_inv + self._decay )                     # dyz orb. of sample
            t += self._coes[..., 6] * self._N_D2 * r_a[..., 0] * (6*r_a_z_sq*r_a_norm_inv**3 + self._decay*(3*r_a_z_sq*r_a_norm_inv_sq-1))                               # dz2 orb. of sample
            t += self._coes[..., 7] * self._N_D * r_a[..., 2] * r_a_norm_inv*( 2*r_a_x_sq*r_a_norm_inv_sq + self._decay*r_a_x_sq*r_a_norm_inv - 1 )                      # dxz orb. of sample
            t += self._coes[..., 8] * self._N_D * r_a_norm_inv*r_a[..., 0]*((r_a_x_sq-r_a_y_sq)*r_a_norm_inv_sq + 0.5*self._decay*(r_a_x_sq-r_a_y_sq)*r_a_norm_inv - 1)  # dx2-y2 orb. of sample
        return t * r_a_norm_inv

    def _missing_conductance(self, tip_orbital: str, *args, **kwargs):
        """Placeholder returning zeros for unsupported tip-sample orbital combinations.

        Raise:
            NotImplementedError
        """
        raise NotImplementedError(f"Not implemented formula for {tip_orbital} tip orb"
                                  f" and {'spd' if self._SAMPLE_ORBITAL_COUNT == 9 else 'sp'} sample orbitals!")

    @property
    @abstractmethod
    def backend(self) -> str:
        """Backend name."""
        pass

    @abstractmethod
    def _to_float(self, array):
        """Cast array to single-precision floating-point number type.

        Parameters
        ----------
        array : array_like
            Input array

        Returns
        -------
        array_like
            a new array of the same shape as the input array, with single-precision floating-point number type.
        """
        pass

    @abstractmethod
    def _unsqueeze(self, array, axis: int | Sequence[int]):
        """Expand the shape of an array by inserting new singleton axes.

        Insert one or more new axes (size 1) at the positions specified by `axis`.

        Parameters
        ----------
        array : array_like
            Input array
        axis : int or sequence of ints
            Position(s) in the expanded shape where new axis/axes are inserted.

        Returns
        -------
        array_like
            View of `array` with the number of dimensions increased by the number of inserted axes
            each inserted dimension has size 1.

        Raises
        ------
        TypeError
            If `axis` has an invalid type.
        ValueError
            If any axis index is out of range for the expanded shape.

        Examples
        --------
        If array.shape == (n_z, n_y, n_x, 3):
        _unsqueeze(array, axis=3) -> (n_z, n_y, n_x, 1, 3)
        _unsqueeze(array, axis=(3, 4)) -> (n_z, n_y, n_x, 1, 1, 3)
        """
        pass

    @abstractmethod
    def _float_zeros(self, size: Sequence[int]):
        """Return a matrix of given shape, filled with single-precision floating-point zeros.

        Parameters
        ----------
        shape : int or sequence of ints
            Shape of the matrix

        Returns
        -------
        matrix
            Single-precision floating-point zero matrix of given shape.
        """
        pass

    @abstractmethod
    def _euclidean_norm(self, array, axis: int):
        """Compute the Euclidean norm along one array axis.

        Parameters
        ----------
        x : array_like
            Input array.
        axis : int
            Axis along which to compute the Euclidean norm.

        Returns
        -------
        float or array_like
            Norm values with `axis` removed from the result shape.
        """
        pass

    @abstractmethod
    def _exp(self, array):
        """Calculate the element-wise exponential of all elements in the input array.

        Parameters
        ----------
        array : array_like
            Input values.

        Returns
        -------
        array_like or scalar
            Output array, element-wise exponential of x. This is a scalar if x is a scalar.
        """
        pass

    @abstractmethod
    def _broadcasted_dot_last_axis(self, a, b):
        """Sum product over the last axis of a and b and broadcast over the uncontracted dimensions.

        Parameters
        ----------
        a : array_like
            First argument.

        b : array_like
            Second argument.

        Returns
        -------
        array_like
            Returns the sum product over the last axis of a and b.
            If a and b are both scalars or both 1-D arrays then a scalar is returned; otherwise an array is returned.
        """
        pass

    @abstractmethod
    def _vstack(self, arrays: List | Tuple):
        """Stack arrays in sequence vertically (row wise).

        This is equivalent to concatenation along the first axis
        after all 1-D arrays of shape `(N,)` have been reshaped to `(1,N)`.

        Parameters
        ----------
        arrays : sequence of array_like
            sequence of arrays to concatenate

        Returns
        -------
        array_like
            The array formed by stacking the given arrays, will be at least 2-D.

        Examples
        --------
        >>> a = [1, 2, 3]
        >>> b = [4, 5, 6]
        >>> self._vstack((a,b))
        tensor([[1, 2, 3],
                [4, 5, 6]])
        >>> a = [[1],[2],[3]]
        >>> b = [[4],[5],[6]]
        >>> self._vstack((a,b))
        tensor([[1],
                [2],
                [3],
                [4],
                [5],
                [6]])
        """
        pass

    @staticmethod
    def _lorentzian(x, loc: float, scale: float):
        """Probability density function (PDF) of the Lorentz or Cauchy distribution.

        Args:
            x (array-like): Evaluation points.
            loc (float): Center of the peak.
            scale (float): Half-width at half-maximum (HWHM), half interquartile range, probable error. Must be > 0.

        Returns:
            array-like: Lorentzian density, same shape as `x`.
        """
        return scale / (math.pi * ((x - loc)**2 + scale ** 2))


class ProbeStmVectorizedSp(ProbeStmVectorized, ABC):
    _SAMPLE_ORBITAL_COUNT = 4

    @property
    def g_tip_orb_fns(self) -> List[Callable]:
        """Return a list of conductance functions for sp sample orbitals,
        where each function represents a different tip orbital.

        Returns:
            List[Callable]: functions computing orbital conductances.
                Each function must return arrays of shape (n_c, n_e, n_a).
        """
        return [
            self._s, self._py, self._pz, self._px,
            partial(self._missing_conductance, tip_orbital="dxy"),
            self._dyz, self._dz2, self._dxz,
            partial(self._missing_conductance, tip_orbital="dx2-y2"),
        ]

    def _dxz(self, from_i, to_i):
        """Orbital conductances for different sp orbitals of sample on a dxz orbital tip.

        Returns:
            array_like, shape (n_c, n_e, n_a)
        """
        r_a = self._r_a[from_i:to_i]  # shape (n_c, 1, n_a, 3)
        r_a_norm = self._r_a_norm[from_i:to_i]  # shape (n_c, 1, n_a)
        decay_sq = self._decay**2
        r_a_norm_inv = 1 / r_a_norm
        r_a_norm_inv_sq = r_a_norm_inv**2
        r_a_norm_inv_cu = r_a_norm_inv**3
        r_a_x_sq = r_a[..., 0]**2
        r_a_z_sq = r_a[..., 2]**2
        r_a_x_mul_z = r_a[..., 0]*r_a[..., 2]
        t =  self._coes[..., 0] * r_a_x_mul_z*(r_a_norm_inv + self._decay)*self._decay  # s  orb. of sample, shape (n_c, n_e, n_a)
        t += self._coes[..., 1] * self._N_P*r_a[..., 1]*r_a_x_mul_z*r_a_norm_inv*( 3*r_a_norm_inv_sq + 3*self._decay*r_a_norm_inv + decay_sq )                                                 # py orb. of sample
        t += self._coes[..., 2] * self._N_P*r_a[..., 0]*( 3*r_a_z_sq*r_a_norm_inv_cu + 3*self._decay*r_a_z_sq*r_a_norm_inv_sq - r_a_norm_inv + decay_sq*r_a_z_sq*r_a_norm_inv - self._decay )  # pz orb. of sample
        t += self._coes[..., 3] * self._N_P*r_a[..., 2]*( 3*r_a_x_sq*r_a_norm_inv_cu + 3*self._decay*r_a_x_sq*r_a_norm_inv_sq - r_a_norm_inv + decay_sq*r_a_x_sq*r_a_norm_inv - self._decay )  # px orb. of sample
        return t * r_a_norm_inv_sq

    def _dyz(self, from_i, to_i):
        """Orbital conductances for different sp orbitals of sample on a dyz orbital tip.

        Returns:
            array_like, shape (n_c, n_e, n_a)
        """
        r_a = self._r_a[from_i:to_i]  # shape (n_c, 1, n_a, 3)
        r_a_norm = self._r_a_norm[from_i:to_i]  # shape (n_c, 1, n_a)
        decay_sq = self._decay**2
        r_a_norm_inv = 1 / r_a_norm
        r_a_norm_inv_sq = r_a_norm_inv**2
        r_a_norm_inv_cu = r_a_norm_inv**3
        r_a_y_sq = r_a[..., 1]**2
        r_a_z_sq = r_a[..., 2]**2
        r_a_y_mul_z = r_a[..., 1]*r_a[..., 2]
        t =  self._coes[..., 0] * r_a_y_mul_z*(r_a_norm_inv + self._decay)*self._decay  # s  orb. of sample, shape (n_c, n_e, n_a)
        t += self._coes[..., 1] * self._N_P*r_a[..., 2]*( 3*r_a_y_sq*r_a_norm_inv_cu + 3*self._decay*r_a_y_sq*r_a_norm_inv_sq - r_a_norm_inv + decay_sq*r_a_y_sq*r_a_norm_inv - self._decay )  # py orb. of sample
        t += self._coes[..., 2] * self._N_P*r_a[..., 1]*( 3*r_a_z_sq*r_a_norm_inv_cu + 3*self._decay*r_a_z_sq*r_a_norm_inv_sq - r_a_norm_inv + decay_sq*r_a_z_sq*r_a_norm_inv - self._decay )  # pz orb. of sample
        t += self._coes[..., 3] * self._N_P*r_a[..., 0]*r_a_y_mul_z*r_a_norm_inv*( 3*r_a_norm_inv_sq + 3*self._decay*r_a_norm_inv + decay_sq )                                                 # px orb. of sample
        return t * r_a_norm_inv_sq

    def _dz2(self, from_i, to_i):
        """Orbital conductances for different sp orbitals of sample on a dz2 orbital tip.

        Returns:
            array_like, shape (n_c, n_e, n_a)
        """
        r_a = self._r_a[from_i:to_i]  # shape (n_c, 1, n_a, 3)
        r_a_norm = self._r_a_norm[from_i:to_i]  # shape (n_c, 1, n_a)
        n_p_div_r_a_norm = self._N_P / r_a_norm
        decay_div_r_a_norm = self._decay / r_a_norm
        decay_sq = self._decay**2
        r_a_norm_inv = 1 / r_a_norm
        r_a_norm_inv_sq = r_a_norm_inv**2
        r_a_norm_inv_cu = r_a_norm_inv**3
        r_a_norm_inv_qu = r_a_norm_inv**4
        r_a_z_sq = r_a[..., 2]**2
        t =  self._coes[..., 0] * (-self._I_3 * decay_sq + self._decay * r_a_z_sq * r_a_norm_inv_cu + decay_sq * r_a_z_sq * r_a_norm_inv_sq - decay_div_r_a_norm)  # s  orb. of sample, shape (n_c, n_e, n_a)
        t += self._coes[..., 1] * n_p_div_r_a_norm*r_a[..., 1]*( 3*r_a_z_sq*r_a_norm_inv_qu + 3*self._decay*r_a_z_sq*r_a_norm_inv_cu - r_a_norm_inv_sq + decay_sq*r_a_z_sq*r_a_norm_inv_sq - decay_div_r_a_norm - self._I_3*decay_sq)      # py orb. of sample
        t += self._coes[..., 2] * n_p_div_r_a_norm*r_a[..., 2]*( 3*r_a_z_sq*r_a_norm_inv_qu + 3*self._decay*r_a_z_sq*r_a_norm_inv_cu - 3*r_a_norm_inv_sq + decay_sq*r_a_z_sq*r_a_norm_inv_sq - 3*decay_div_r_a_norm - self._I_3*decay_sq)  # pz orb. of sample
        t += self._coes[..., 3] * n_p_div_r_a_norm*r_a[..., 0]*( 3*r_a_z_sq*r_a_norm_inv_qu + 3*self._decay*r_a_z_sq*r_a_norm_inv_cu - r_a_norm_inv_sq + decay_sq*r_a_z_sq*r_a_norm_inv_sq - decay_div_r_a_norm - self._I_3*decay_sq)      # px orb. of sample
        return t


class ProbeStmVectorizedSpd(ProbeStmVectorized, ABC):
    _SAMPLE_ORBITAL_COUNT = 9

    @property
    def g_tip_orb_fns(self) -> List[Callable]:
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
