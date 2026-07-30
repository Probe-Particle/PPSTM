from abc import ABC
from typing import Callable

import numpy as np


class _TestDidv(ABC):
    """Base test class for dI/dV (differential conductance) calculations.

    Validates the correctness of dI/dV backend implementations for STM simulations.
    """
    _SAMPLE_ORBITAL_COUNT: int
    _V: float  # voltage
    _EIGENENERGIES: np.ndarray
    _LCAO_COEFFICIENTS: np.ndarray
    _SAMPLE_ATOMS_POSITIONS: np.ndarray

    _WORK_FUNCTION = 5.  # eV
    _ETA = 1e-7  # eV, Lorentzian width of states in energy scale

    _SEED = 3

    def test_n_eigenenergies_eq_n_lcao_rows(self):
        """Check that the size of the 1st dimension of the LCAO coefficients equals the nr. eigenenergies."""
        assert len(self._LCAO_COEFFICIENTS) == len(self._EIGENENERGIES)

    def test_n_lcao_cols_eq_n_sample_positions_times_n_sample_orbitals(self):
        """Check that the size of the 2nd dimension of the LCAO coefficients equals
        the nr. sample atoms times the nr. sample orbitals."""
        assert self._LCAO_COEFFICIENTS.shape[1] == len(self._SAMPLE_ATOMS_POSITIONS) * self._SAMPLE_ORBITAL_COUNT

    def _didv_generic(self, tip_positions: np.ndarray, tip_orb: np.ndarray, backend: Callable) -> np.ndarray:
        """Execute dI/dV calculation with the given backend implementation.

        Args:
            tip_positions (np.ndarray): array of tip coordinates, shape: (n_z, n_y, n_x, 3)
            tip_orb (np.ndarray): tip-orbital coefficients, shape: (9,)
            backend (Callable): dI/dV implementation

        Returns:
            dI/dV signal array from the backend, shape: (n_z, n_y, n_x)
        """
        return backend(V=self._V,
                       WF=self._WORK_FUNCTION,
                       eta=self._ETA,
                       eig=self._EIGENENERGIES,
                       R=tip_positions,
                       Rat=self._SAMPLE_ATOMS_POSITIONS,
                       coes=self._LCAO_COEFFICIENTS,
                       tip_coes=tip_orb,
                       orb_t=self._SAMPLE_ORBITAL_COUNT)

    def _make_rng(self) -> np.random.Generator:
        """Return a deterministic random number generator for reproducible test data."""
        return np.random.default_rng(self._SEED)
