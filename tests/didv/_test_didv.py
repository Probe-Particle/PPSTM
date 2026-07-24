from abc import ABC
from typing import Callable

import numpy as np


class _TestDidv(ABC):
    """Base test class for dIdV (differential conductance) calculations.

    Validates the correctness of dIdV backend implementations for STM simulations.
    """
    _SAMPLE_ORBITAL_COUNT: int
    _V: float  # voltage
    _ETA: float  # eV, Lorentzian width of states in energy scale
    _EIGENENERGIES: np.ndarray
    _LCAO_KS_COEFFICIENTS: np.ndarray
    _SAMPLE_ATOMS_POSITIONS: np.ndarray

    _WORK_FUNCTION = 5.  # eV

    _SEED = 3

    def test_n_eigenenergies_eq_n_lcao_rows(self):
        assert len(self._EIGENENERGIES) == len(self._LCAO_KS_COEFFICIENTS)

    def test_n_lcao_cols_eq_n_sample_positions_times_n_sample_orbitals(self):
        assert self._LCAO_KS_COEFFICIENTS.shape[1] == self._SAMPLE_ORBITAL_COUNT * len(self._SAMPLE_ATOMS_POSITIONS)

    def _didv_generic(self, tip_positions: np.ndarray, tip_orb: np.ndarray, backend: Callable) -> np.ndarray:
        """Execute dIdV calculation with the given backend implementation.

        Args:
            tip_positions (np.ndarray): array of tip coordinates, shape: (n_z, n_y, n_x, 3)
            tip_orb (np.ndarray): conductance coefficient by tip orbital, shape: (9,)
            backend (Callable): dI/dV implementation

        Returns:
            dIdV signal array from the backend, shape: (n_z, n_y, n_x)
        """
        return backend(V=self._V,
                       WF=self._WORK_FUNCTION,
                       eta=self._ETA,
                       eig=self._EIGENENERGIES,
                       R=tip_positions,
                       Rat=self._SAMPLE_ATOMS_POSITIONS,
                       coes=self._LCAO_KS_COEFFICIENTS,
                       tip_coes=tip_orb,
                       orb_t=self._SAMPLE_ORBITAL_COUNT)

    def _make_rng(self) -> np.random.Generator:
        return np.random.default_rng(self._SEED)
