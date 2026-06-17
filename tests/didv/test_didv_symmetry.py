"""Inspired by examples/orbitals_tests/orbitals.toml."""
import itertools
from abc import ABC
from typing import Callable

import numpy as np
import pytest

from pyPPSTM.probe_stm_opencl import ProbeSTMOpenCLParallel, ProbeSTMOpenCLSequential
from pyPPSTM.probe_stm_numpy import ProbeStmNumpy
from pyPPSTM.probe_stm_pytorch import ProbeStmPytorch
from tests.didv import _test_didv


class _TestDidvSymmetry(_test_didv._TestDidv, ABC):
    _LOWEST_TIP_POSITION = -5.0
    _HIGHEST_TIP_POSITION = 5.0
    _TIP_POSITION_STEP = 0.2
    _NR_TIP_POSITIONS_SINGLE_AXIS = int((_HIGHEST_TIP_POSITION - _LOWEST_TIP_POSITION) / _TIP_POSITION_STEP) + 1
    _POSITIVE_TIP_POSITIONS_SINGLE_AXIS = (
        (lambda t: t[t>0.0])
        (np.linspace(start=_LOWEST_TIP_POSITION,
                     stop=_HIGHEST_TIP_POSITION,
                     num=_NR_TIP_POSITIONS_SINGLE_AXIS))
    )
    _NR_POSITIVE_TIP_POSITIONS_SINGLE_AXIS = len(_POSITIVE_TIP_POSITIONS_SINGLE_AXIS)
    _NR_SAMPLED_TIP_POSITIONS = 10 # tot: 25*25*25

    _V = -0.04   # voltage
    _ETA = 1e-7  # eV, Lorentzian width of states in energy scale
    _TIP = np.asarray([1., ] + [0.0, ] * 8) # TIP_ORB = 's'

    _SAMPLE_ATOMS_POSITIONS = np.asarray([[0.0, 0.0, 0.0]])

    def test_even_symmetry_at_tip_center_for_single_tip_position(self, tip_position: np.ndarray,
                                                                 didv_backend: Callable,
                                                                 rtol:float):
        """Verify dIdV curve is even-symmetric about the tip's center position for a single tip positions.
        Checks that f(x) == f(-x) for all sampled tip positions."""
        self._test_even_symmetry_at_tip_center(tip_position, didv_backend, rtol)

    def test_even_symmetry_at_tip_center_for_tip_positions_batch(self,
                                                                 all_tip_positions_batch: np.ndarray,
                                                                 didv_backend: Callable,
                                                                 rtol:float):
        """Verify dIdV curve is even-symmetric about the tip's center position for a batch of tip positions.
        Checks that f(x) == f(-x) for all sampled tip positions."""
        self._test_even_symmetry_at_tip_center(all_tip_positions_batch, didv_backend, rtol)

    def _test_even_symmetry_at_tip_center(self,
                                          tip_positions: np.ndarray,
                                          didv_backend: Callable,
                                          rtol: float):
        """Verify dIdV curve is even-symmetric about the tip's center position.
        Checks that f(x) == f(-x) for all sampled tip positions."""
        didv = self._didv(tip_positions, backend=didv_backend)
        didv_with_opposite_tip_position = self._didv(-tip_positions, backend=didv_backend)

        np.testing.assert_allclose(didv, didv_with_opposite_tip_position, rtol=rtol)

    def _didv(self, tip_positions: np.ndarray, backend: Callable) -> np.ndarray:
        return super()._didv_generic(tip_positions=tip_positions,
                                     tip_orb=self._TIP,
                                     backend=backend)

    @pytest.fixture(scope="class")
    def all_tip_positions_batch(self):
        yield (np.asarray(list(itertools.product(*(self._POSITIVE_TIP_POSITIONS_SINGLE_AXIS,) * 3)))
               * self._make_rng().choice([-1.0, 1.0], size=(self._NR_POSITIVE_TIP_POSITIONS_SINGLE_AXIS**3, 3)))\
            .reshape(*(self._NR_POSITIVE_TIP_POSITIONS_SINGLE_AXIS,)*3, 3)

    def pytest_generate_tests(self, metafunc):
        if "tip_position" in metafunc.fixturenames:
            tip_positions = list(self._sample_tip_positions(self._NR_SAMPLED_TIP_POSITIONS)
                                 .reshape(self._NR_SAMPLED_TIP_POSITIONS, 1, 1, 1, 3))
            metafunc.parametrize("tip_position", tip_positions,
                                 ids=map(lambda x: str(x.reshape(3)), tip_positions))

        if "didv_backend" in metafunc.fixturenames and \
           "rtol" in metafunc.fixturenames:
            metafunc.parametrize(("didv_backend", "rtol"),
                                 (
                                         (ProbeSTMOpenCLParallel.didv,   1e-7),
                                         (ProbeSTMOpenCLSequential.didv, 1e-7),
                                         (ProbeStmNumpy.didv,                           1e-7),
                                         (ProbeStmPytorch.didv,                         1e-7),
                                 ),
                                 ids=[
                                     "OpenCL parallel-rtol=1e-7",
                                     "OpenCL sequential-rtol=1e-7",
                                     "NumPy-rtol=1e-7",
                                     "PyTorch-rtol=1e-7",
                                 ])

    def _sample_tip_positions(self, n_samples: int) -> np.ndarray:
        rng = self._make_rng()
        positive_tip_positions_indices = rng.choice(range(self._NR_POSITIVE_TIP_POSITIONS_SINGLE_AXIS**3),
                                                    size=n_samples, replace=False)
        x_i, y_i, z_i = np.unravel_index(positive_tip_positions_indices,
                                         shape=(self._NR_POSITIVE_TIP_POSITIONS_SINGLE_AXIS,)*3)
        positive_tip_positions = np.stack((
            self._POSITIVE_TIP_POSITIONS_SINGLE_AXIS[x_i],
            self._POSITIVE_TIP_POSITIONS_SINGLE_AXIS[y_i],
            self._POSITIVE_TIP_POSITIONS_SINGLE_AXIS[z_i],
        ), axis=-1)
        sign_tip_positions = rng.choice([-1.0, 1.0], size=(n_samples, 3))
        return positive_tip_positions * sign_tip_positions


class TestDidvSymmetrySp(_TestDidvSymmetry):
    _SAMPLE_ORBITAL_COUNT = 4
    _EIGENENERGIES = np.asarray([-0.04, -0.03, -0.02, -0.01])
    _LCAO_KS_COEFFICIENTS_DIAG = [1.0, 1.0, 1.0, 1.0]
    _LCAO_KS_COEFFICIENTS = np.diag(_LCAO_KS_COEFFICIENTS_DIAG)


class TestDidvSymmetrySpd(_TestDidvSymmetry):
    _SAMPLE_ORBITAL_COUNT = 9
    _EIGENENERGIES = np.asarray([-0.04, -0.03, -0.02, -0.01, 0.0, 0.01, 0.02, 0.03, 0.04])
    _LCAO_KS_COEFFICIENTS_DIAG = [1.0, 1.0, 1.0, 1.0, 0.2, 0.2, 0.2, 0.2, 0.2]
    _LCAO_KS_COEFFICIENTS = np.diag(_LCAO_KS_COEFFICIENTS_DIAG)
