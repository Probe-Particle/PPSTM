"""
STM (Scanning Tunneling Microscopy) calculations using vectorization.

This module provides vectorized implementations for computing dI/dV (differential
conductance) spectra with sp(d)-sp(d) orbital interactions.
"""
import logging
from abc import abstractmethod, ABC
from functools import partial
from typing import List, Callable, Sequence

import numpy as np
import math


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
                 tip_coes):
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
        """
        v = np.float32(V)
        wf = np.float32(WF)
        eta = np.float32(eta)
        self._tip_dos = self._to_float(tip_coes)  # shape (9,)
        eig = self._to_float(eig)                                                     # shape                (n_e,)
        self._sample_dos = self._lorentzian(x=eig, loc=v, scale=0.5 * eta)            # shape                (n_e,)

        sample_atom_position = self._to_float(Rat)                                    # shape                     (n_a, 3)
        tip_position = self._unsqueeze(self._to_float(R), axis=[3, 4])                # shape (n_z, n_y, n_x,   1,   1, 3)

        self._decay = math.sqrt((2 * wf) * self._EV)
        self._r_a = (tip_position - sample_atom_position) * self._AB                  # shape (n_z, n_y, n_x,   1, n_a, 3)
        self._r_a_norm = self._norm(self._r_a, ord=2, axis=-1)                        # shape (n_z, n_y, n_x,   1, n_a)
        self._radial = self._exp(-self._decay * self._r_a_norm)                       # shape (n_z, n_y, n_x,   1, n_a)

        self._coes = self._to_float(coes)\
                         .reshape(len(eig),
                                  len(Rat),
                                  self._SAMPLE_ORBITAL_COUNT)                         # shape                (n_e, n_a,    orb_t)

        self._logger = logging.getLogger(self.__class__.__name__)

    def __call__(self):
        """Compute dI/dV."""

        self._logger.debug(f"Entering the dI/dV ( sp(d)-sp(d) ) {self.backend} procedure")

        g_tip_orbs_sum = self._float_zeros((*self._r_a.shape[:3], self._coes.shape[0]))
        for i, g_tip_orb_fn in enumerate(self.g_tip_orb_fns):
            if self._tip_dos[i] > 0.:
                g_tip_orb = self._einsum("...a,...a",
                                         g_tip_orb_fn(),      # shape (n_z, n_y, n_x, n_e, n_a)
                                         self._radial)                   # shape (n_z, n_y, n_x, n_e)

                g_tip_orbs_sum += self._tip_dos[i] * g_tip_orb ** 2

        g = self._einsum("...i,i",
                         g_tip_orbs_sum,
                         self._sample_dos)                                                    # shape (n_z, n_y, n_x)

        g *= 16 * math.pi ** 3 * self._decay

        self._logger.debug(f"dI/dV ( sp(d)-sp(d) ) {self.backend} procedure DONE")
        return g

    @property
    @abstractmethod
    def g_tip_orb_fns(self) -> List[Callable]:
        """Return a list of conductance functions for specific sample orbitals,
        where each function represents a different tip orbital.

        Returns:
            List[Callable]: functions computing orbital conductances.
                Each function must return arrays of shape (n_z, n_y, n_x, n_e, n_a).
        """
        pass

    def _s(self):
        """Orbital conductances for different sp or spd orbitals of sample on an s orbital tip.

        Returns:
            array_like, shape (n_z, n_y, n_x, n_e, n_a)
        """
        n_p_div_r_a_norm = self._N_P / self._r_a_norm  # shape (n_z, n_y, n_x, 1, n_a)
        r_a_norm_sq = self._r_a_norm ** 2
        n_d_div_r_a_norm_sq = self._N_D / r_a_norm_sq

        t = self._coes[..., 0]                                             # s  orb. of sample, shape                (n_e, n_a)
        t = t + self._coes[..., 1] * n_p_div_r_a_norm * self._r_a[..., 1]  # py orb. of sample, shape (n_z, n_y, n_x, n_e, n_a)
        t = t + self._coes[..., 2] * n_p_div_r_a_norm * self._r_a[..., 2]  # pz orb. of sample
        t = t + self._coes[..., 3] * n_p_div_r_a_norm * self._r_a[..., 0]  # px orb. of sample
        if self._SAMPLE_ORBITAL_COUNT == 9:
            t = t + self._coes[..., 4] * n_d_div_r_a_norm_sq * self._r_a[..., 0] * self._r_a[..., 1]  # dxy orb. of sample
            t = t + self._coes[..., 5] * n_d_div_r_a_norm_sq * self._r_a[..., 1] * self._r_a[..., 2]  # dyz orb. of sample
            t = t + self._coes[..., 6] * self._N_D2 * (3 * self._r_a[..., 2] ** 2 / r_a_norm_sq - 1)  # dz2 orb. of sample
            t = t + self._coes[..., 7] * n_d_div_r_a_norm_sq * self._r_a[..., 0] * self._r_a[..., 2]  # dxz orb. of sample
            t = t + self._coes[..., 8] * n_d_div_r_a_norm_sq * .5 * (self._r_a[..., 0] ** 2 - self._r_a[..., 1] ** 2)  # dx2-y2 orb. of sample
        return t

    def _py(self):
        """Orbital conductances for different sp or spd orbitals of sample on a py orbital tip.

        Returns:
            array_like, shape (n_z, n_y, n_x, n_e, n_a)
        """
        r_a_norm_inv = 1 / self._r_a_norm
        r_a_norm_inv_sq = r_a_norm_inv**2
        r_a_x_sq = self._r_a[..., 0]**2
        r_a_y_sq = self._r_a[..., 1]**2
        r_a_z_sq = self._r_a[..., 2]**2
        t =  self._coes[..., 0] * self._r_a[..., 1] * self._decay  # s  orb. of sample, shape (n_z, n_y, n_x, n_e, n_a)
        t += self._coes[..., 1] * self._N_P * ( -1 + self._decay*r_a_norm_inv*r_a_y_sq + r_a_norm_inv_sq*r_a_y_sq )     # py orb. of sample
        t += self._coes[..., 2] * self._N_P * self._r_a[..., 1] * self._r_a[..., 2] * (self._decay*r_a_norm_inv + r_a_norm_inv_sq)  # pz orb. of sample
        t += self._coes[..., 3] * self._N_P * self._r_a[..., 1] * self._r_a[..., 0] * (self._decay*r_a_norm_inv + r_a_norm_inv_sq)  # px orb. of sample
        if self._SAMPLE_ORBITAL_COUNT == 9:
            t += self._coes[..., 4] * self._N_D * self._r_a[..., 0] * r_a_norm_inv*( 2*r_a_y_sq*r_a_norm_inv_sq + self._decay*r_a_y_sq*r_a_norm_inv - 1 )   # dxy orb. of sample
            t += self._coes[..., 5] * self._N_D * self._r_a[..., 2] * r_a_norm_inv*( 2*r_a_y_sq*r_a_norm_inv_sq + self._decay*r_a_y_sq*r_a_norm_inv - 1 )   # dyz orb. of sample
            t += self._coes[..., 6] * self._N_D2 * self._r_a[..., 1] * (6*r_a_z_sq*r_a_norm_inv**3 + self._decay*(3*r_a_z_sq*r_a_norm_inv_sq-1))            # dz2 orb. of sample
            t += self._coes[..., 7] * self._N_D * self._r_a[..., 0] * self._r_a[..., 1]*self._r_a[..., 2]*r_a_norm_inv_sq*( 2*r_a_norm_inv + self._decay )  # dxz orb. of sample
            t += self._coes[..., 8] * self._N_D * r_a_norm_inv*self._r_a[..., 1]*((r_a_x_sq-r_a_y_sq)*r_a_norm_inv_sq + 0.5*self._decay*(r_a_x_sq-r_a_y_sq)*r_a_norm_inv + 1)  # dx2-y2 orb. of sample
        return t * r_a_norm_inv

    def _pz(self):
        """Orbital conductances for different sp or spd orbitals of sample on a pz orbital tip.

        Returns:
            array_like, shape (n_z, n_y, n_x, n_e, n_a)
        """
        r_a_norm_inv = 1 / self._r_a_norm
        r_a_norm_inv_sq = r_a_norm_inv**2
        r_a_z_sq = self._r_a[..., 2]**2
        t =  self._coes[..., 0] * self._r_a[..., 2] * self._decay  # s  orb. of sample, shape (n_z, n_y, n_x, n_e, n_a)
        t += self._coes[..., 1] * self._N_P * self._r_a[..., 2] * self._r_a[..., 1] * (self._decay*r_a_norm_inv + r_a_norm_inv_sq)  # py orb. of sample
        t += self._coes[..., 2] * self._N_P * ( -1 + self._decay*r_a_norm_inv*r_a_z_sq + r_a_norm_inv_sq*r_a_z_sq )                 # pz orb. of sample
        t += self._coes[..., 3] * self._N_P * self._r_a[..., 2] * self._r_a[..., 0]*( self._decay*r_a_norm_inv + r_a_norm_inv_sq )  # px orb. of sample
        if self._SAMPLE_ORBITAL_COUNT == 9:
            t += self._coes[..., 4] * self._N_D * self._r_a[..., 0] * self._r_a[..., 1]*self._r_a[..., 2]*r_a_norm_inv_sq*( 2*r_a_norm_inv + self._decay )                   # dxy orb. of sample
            t += self._coes[..., 5] * self._N_D * self._r_a[..., 1] * r_a_norm_inv*( 2*r_a_z_sq*r_a_norm_inv_sq + self._decay*r_a_z_sq*r_a_norm_inv - 1 )                    # dyz orb. of sample
            t += self._coes[..., 6] * self._N_D2 * ((6*self._r_a[..., 2]**3*r_a_norm_inv**3-6*self._r_a[..., 2]*r_a_norm_inv) + self._decay*(3*r_a_z_sq*r_a_norm_inv_sq-1))  # dz2 orb. of sample
            t += self._coes[..., 7] * self._N_D * self._r_a[..., 0] * r_a_norm_inv*( 2*r_a_z_sq*r_a_norm_inv_sq + self._decay*r_a_z_sq*r_a_norm_inv - 1 )                    # dxz orb. of sample
            t += self._coes[..., 8] * self._N_D * r_a_norm_inv_sq*self._r_a[..., 2]*(self._r_a[..., 0]**2-self._r_a[..., 1]**2)*( r_a_norm_inv + 0.5*self._decay  )          # dx2-y2 orb. of sample
        return t * r_a_norm_inv

    def _px(self):
        """Orbital conductances for different sp or spd orbitals of sample on a px orbital tip.

        Returns:
            array_like, shape (n_z, n_y, n_x, n_e, n_a)
        """
        r_a_norm_inv = 1 / self._r_a_norm
        r_a_norm_inv_sq = r_a_norm_inv**2
        r_a_x_sq = self._r_a[..., 0]**2
        r_a_y_sq = self._r_a[..., 1]**2
        r_a_z_sq = self._r_a[..., 2]**2
        t =  self._coes[..., 0] * self._r_a[..., 0] * self._decay  # s  orb. of sample, shape (n_z, n_y, n_x, n_e, n_a)
        t += self._coes[..., 1] * self._N_P * self._r_a[..., 0] * self._r_a[..., 1] * (self._decay*r_a_norm_inv + r_a_norm_inv_sq)  # py orb. of sample
        t += self._coes[..., 2] * self._N_P * self._r_a[..., 0] * self._r_a[..., 2] * (self._decay*r_a_norm_inv + r_a_norm_inv_sq)  # pz orb. of sample
        t += self._coes[..., 3] * self._N_P * ( -1 + self._decay*r_a_norm_inv*r_a_x_sq + r_a_norm_inv_sq*r_a_x_sq )                 # px orb. of sample
        if self._SAMPLE_ORBITAL_COUNT == 9:
            t += self._coes[..., 4] * self._N_D * self._r_a[..., 1] * r_a_norm_inv*( 2*r_a_x_sq*r_a_norm_inv_sq + self._decay*r_a_x_sq*r_a_norm_inv - 1 )                      # dxy orb. of sample
            t += self._coes[..., 5] * self._N_D * self._r_a[..., 0] * self._r_a[..., 1]*self._r_a[..., 2]*r_a_norm_inv_sq*( 2*r_a_norm_inv + self._decay )                     # dyz orb. of sample
            t += self._coes[..., 6] * self._N_D2 * self._r_a[..., 0] * (6*r_a_z_sq*r_a_norm_inv**3 + self._decay*(3*r_a_z_sq*r_a_norm_inv_sq-1))                               # dz2 orb. of sample
            t += self._coes[..., 7] * self._N_D * self._r_a[..., 2] * r_a_norm_inv*( 2*r_a_x_sq*r_a_norm_inv_sq + self._decay*r_a_x_sq*r_a_norm_inv - 1 )                      # dxz orb. of sample
            t += self._coes[..., 8] * self._N_D * r_a_norm_inv*self._r_a[..., 0]*((r_a_x_sq-r_a_y_sq)*r_a_norm_inv_sq + 0.5*self._decay*(r_a_x_sq-r_a_y_sq)*r_a_norm_inv - 1)  # dx2-y2 orb. of sample
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
    def _norm(self, array, ord: int, axis: int):
        """Matrix or vector norm.

        This function is able to return one of eight different matrix norms,
        or one of an infinite number of vector norms (described below), depending
        on the value of the ``ord`` parameter.

        Parameters
        ----------
        x : array_like
            Input array.  If `axis` is None, `x` must be 1-D or 2-D, unless `ord`
            is None. If both `axis` and `ord` are None, the 2-norm of
            ``x.ravel`` will be returned.
        ord : {int, float, inf, -inf, 'fro', 'nuc'}, optional
            Order of the norm (see table under ``Notes`` for what values are
            supported for matrices and vectors respectively). inf means numpy's
            `inf` object. The default is None.
        axis : {None, int, 2-tuple of ints}, optional.
            If `axis` is an integer, it specifies the axis of `x` along which to
            compute the vector norms.  If `axis` is a 2-tuple, it specifies the
            axes that hold 2-D matrices, and the matrix norms of these matrices
            are computed.  If `axis` is None then either a vector norm (when `x`
            is 1-D) or a matrix norm (when `x` is 2-D) is returned. The default
            is None.

        Returns
        -------
        float or array_like
            Norm of the matrix or vector(s).

        Notes
        -----
        The following norms can be calculated:

        =====  ============================  ==========================
        ord    norm for matrices             norm for vectors
        =====  ============================  ==========================
        None   Frobenius norm                2-norm
        'fro'  Frobenius norm                --
        'nuc'  nuclear norm                  --
        inf    max(sum(abs(x), axis=1))      max(abs(x))
        -inf   min(sum(abs(x), axis=1))      min(abs(x))
        0      --                            sum(x != 0)
        1      max(sum(abs(x), axis=0))      as below
        -1     min(sum(abs(x), axis=0))      as below
        2      2-norm (largest sing. value)  as below
        -2     smallest singular value       as below
        other  --                            sum(abs(x)**ord)**(1./ord)
        =====  ============================  ==========================
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
    def _einsum(self, subscripts, *operands):
        """Evaluates the Einstein summation convention on the operands.

        Using the Einstein summation convention, many common multi-dimensional,
        linear algebraic array operations can be represented in a simple fashion.
        In *implicit* mode `einsum` computes these values.

        In *explicit* mode, `einsum` provides further flexibility to compute
        other array operations that might not be considered classical Einstein
        summation operations, by disabling, or forcing summation over specified
        subscript labels.

        Parameters
        ----------
        subscripts : str
            Specifies the subscripts for summation as comma separated list of
            subscript labels. An implicit (classical Einstein summation)
            calculation is performed unless the explicit indicator '->' is
            included as well as subscript labels of the precise output form.
        operands : list of array_like
            These are the arrays for the operation.

        Returns
        -------
        array_like
            The calculation based on the Einstein summation convention.
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
                Each function must return arrays of shape (n_z, n_y, n_x, n_e, n_a).
        """
        return [
            self._s, self._py, self._pz, self._px,
            partial(self._missing_conductance, tip_orbital="dxy"),
            self._dyz, self._dz2, self._dxz,
            partial(self._missing_conductance, tip_orbital="dx2-y2"),
        ]

    def _dxz(self):
        """Orbital conductances for different sp orbitals of sample on a dxz orbital tip.

        Returns:
            array_like, shape (n_z, n_y, n_x, n_e, n_a)
        """
        decay_sq = self._decay**2
        r_a_norm_inv = 1 / self._r_a_norm
        r_a_norm_inv_sq = r_a_norm_inv**2
        r_a_norm_inv_cu = r_a_norm_inv**3
        r_a_x_sq = self._r_a[..., 0]**2
        r_a_z_sq = self._r_a[..., 2]**2
        r_a_x_mul_z = self._r_a[..., 0]*self._r_a[..., 2]
        t =  self._coes[..., 0] * r_a_x_mul_z*(r_a_norm_inv + self._decay)*self._decay  # s  orb. of sample, shape (n_z, n_y, n_x, n_e, n_a)
        t += self._coes[..., 1] * self._N_P*self._r_a[..., 1]*r_a_x_mul_z*r_a_norm_inv*( 3*r_a_norm_inv_sq + 3*self._decay*r_a_norm_inv + decay_sq )                                                 # py orb. of sample
        t += self._coes[..., 2] * self._N_P*self._r_a[..., 0]*( 3*r_a_z_sq*r_a_norm_inv_cu + 3*self._decay*r_a_z_sq*r_a_norm_inv_sq - r_a_norm_inv + decay_sq*r_a_z_sq*r_a_norm_inv - self._decay )  # pz orb. of sample
        t += self._coes[..., 3] * self._N_P*self._r_a[..., 2]*( 3*r_a_x_sq*r_a_norm_inv_cu + 3*self._decay*r_a_x_sq*r_a_norm_inv_sq - r_a_norm_inv + decay_sq*r_a_x_sq*r_a_norm_inv - self._decay )  # px orb. of sample
        return t * r_a_norm_inv_sq

    def _dyz(self):
        """Orbital conductances for different sp orbitals of sample on a dyz orbital tip.

        Returns:
            array_like, shape (n_z, n_y, n_x, n_e, n_a)
        """
        decay_sq = self._decay**2
        r_a_norm_inv = 1 / self._r_a_norm
        r_a_norm_inv_sq = r_a_norm_inv**2
        r_a_norm_inv_cu = r_a_norm_inv**3
        r_a_y_sq = self._r_a[..., 1]**2
        r_a_z_sq = self._r_a[..., 2]**2
        r_a_y_mul_z = self._r_a[..., 1]*self._r_a[..., 2]
        t =  self._coes[..., 0] * r_a_y_mul_z*(r_a_norm_inv + self._decay)*self._decay  # s  orb. of sample, shape (n_z, n_y, n_x, n_e, n_a)
        t += self._coes[..., 1] * self._N_P*self._r_a[..., 2]*( 3*r_a_y_sq*r_a_norm_inv_cu + 3*self._decay*r_a_y_sq*r_a_norm_inv_sq - r_a_norm_inv + decay_sq*r_a_y_sq*r_a_norm_inv - self._decay )  # py orb. of sample
        t += self._coes[..., 2] * self._N_P*self._r_a[..., 1]*( 3*r_a_z_sq*r_a_norm_inv_cu + 3*self._decay*r_a_z_sq*r_a_norm_inv_sq - r_a_norm_inv + decay_sq*r_a_z_sq*r_a_norm_inv - self._decay )  # pz orb. of sample
        t += self._coes[..., 3] * self._N_P*self._r_a[..., 0]*r_a_y_mul_z*r_a_norm_inv*( 3*r_a_norm_inv_sq + 3*self._decay*r_a_norm_inv + decay_sq )                                                 # px orb. of sample
        return t * r_a_norm_inv_sq

    def _dz2(self):
        """Orbital conductances for different sp orbitals of sample on a dz2 orbital tip.

        Returns:
            array_like, shape (n_z, n_y, n_x, n_e, n_a)
        """
        n_p_div_r_a_norm = self._N_P / self._r_a_norm
        decay_div_r_a_norm = self._decay / self._r_a_norm
        decay_sq = self._decay**2
        r_a_norm_inv = 1 / self._r_a_norm
        r_a_norm_inv_sq = r_a_norm_inv**2
        r_a_norm_inv_cu = r_a_norm_inv**3
        r_a_norm_inv_qu = r_a_norm_inv**4
        r_a_z_sq = self._r_a[..., 2]**2
        t =  self._coes[..., 0] * (-self._I_3 * decay_sq + self._decay * r_a_z_sq * r_a_norm_inv_cu + decay_sq * r_a_z_sq * r_a_norm_inv_sq - decay_div_r_a_norm)  # s  orb. of sample, shape (n_z, n_y, n_x, n_e, n_a)
        t += self._coes[..., 1] * n_p_div_r_a_norm*self._r_a[..., 1]*( 3*r_a_z_sq*r_a_norm_inv_qu + 3*self._decay*r_a_z_sq*r_a_norm_inv_cu - r_a_norm_inv_sq + decay_sq*r_a_z_sq*r_a_norm_inv_sq - decay_div_r_a_norm - self._I_3*decay_sq)      # py orb. of sample
        t += self._coes[..., 2] * n_p_div_r_a_norm*self._r_a[..., 2]*( 3*r_a_z_sq*r_a_norm_inv_qu + 3*self._decay*r_a_z_sq*r_a_norm_inv_cu - 3*r_a_norm_inv_sq + decay_sq*r_a_z_sq*r_a_norm_inv_sq - 3*decay_div_r_a_norm - self._I_3*decay_sq)  # pz orb. of sample
        t += self._coes[..., 3] * n_p_div_r_a_norm*self._r_a[..., 0]*( 3*r_a_z_sq*r_a_norm_inv_qu + 3*self._decay*r_a_z_sq*r_a_norm_inv_cu - r_a_norm_inv_sq + decay_sq*r_a_z_sq*r_a_norm_inv_sq - decay_div_r_a_norm - self._I_3*decay_sq)      # px orb. of sample
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
