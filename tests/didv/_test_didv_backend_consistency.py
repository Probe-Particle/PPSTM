from abc import ABC, abstractmethod
from typing import Callable, List

import numpy as np
import pytest

from pyPPSTM import ProbeSTM
from pyPPSTM.probe_stm_numpy import ProbeStmNumpy
from pyPPSTM.probe_stm_opencl import ProbeSTMOpenCLParallel, ProbeSTMOpenCLSequential
from pyPPSTM.probe_stm_pytorch import ProbeStmPytorch
from tests.didv import _test_didv


class _TestDidvBackendConsistency(_test_didv._TestDidv, ABC):
    """Verify that differential conductance (dI/dV) from the C++ backend matches alternative backends.

    This abstract base class provides a comprehensive test suite that validates the C++ reference
    implementation of dI/dV calculations against multiple alternative backends (OpenCL, NumPy)
    across various tip and sample orbitals.

    Inspired by examples/orbitals_tests."""
    _TIP_ORBITALS: List[np.ndarray]
    _TIP_ORBITAL_LABELS: List[str]

    _NR_SAMPLED_TIP_POSITIONS = 10 # tot: 51*51*51

    _V = -.04    # voltage
    _ETA = 1e-7  # eV, Lorentzian width of states in energy scale
    _S_TIP_ORB     = np.asarray([ 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    _PZ_TIP_ORB    = np.asarray([ 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    _PXY_TIP_ORB   = np.asarray([ 0.0, 0.5, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0])
    _CO_TIP_ORB    = np.asarray([0.15, 0.5, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0])
    _DZ2_TIP_ORB   = np.asarray([ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0])
    _DXZYZ_TIP_ORB = np.asarray([ 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.5, 0.0])

    _SAMPLE_ATOMS_POSITIONS = np.asarray([[0., 0., 0.]])  # Atom positions, types, charges

    def test_didv_cpp_matches_other_backends_for_single_tip_position(self,
                                                                     tip_position: np.ndarray,
                                                                     tip_orb: np.ndarray,
                                                                     didv_backend: Callable,
                                                                     rtol: float):
        """Verify that dI/dV from the C++ backend matches other backends for a single tip positions."""
        self._test_didv_cpp_matches_other_backends(tip_position,
                                                   tip_orb,
                                                   didv_backend,
                                                   rtol)

    def test_didv_cpp_matches_other_backends_for_tip_positions_batch(self,
                                                                     all_tip_positions_batch: np.ndarray,
                                                                     tip_orb: np.ndarray,
                                                                     didv_backend: Callable,
                                                                     rtol: float):
        """Verify that dI/dV from the C++ backend matches other backends for a batch of tip positions."""
        self._test_didv_cpp_matches_other_backends(all_tip_positions_batch,
                                                   tip_orb,
                                                   didv_backend,
                                                   rtol)

    def _test_didv_cpp_matches_other_backends(self,
                                              tip_positions: np.ndarray,
                                              tip_orb: np.ndarray,
                                              didv_backend: Callable,
                                              rtol: float):
        """
        Verify that dI/dV from the C++ backend matches other backends for a batch of tip positions.

        This test computes the expected differential conductance (dI/dV) using the C++
        reference implementation and compares it to the value
        returned by the parametrized backend.
        Comparison is performed over a batch of tip positions and a single tip orbital
        with the supplied relative tolerance (rtol).

        Args:
            tip_positions: (np.ndarray): batch of tip positions, shape (n_z, n_y, n_x, 3).
            tip_orb (np.ndarray): conductance coefficient by tip orbital, shape (9,).
            didv_backend (callable): dI/dV implementation under test.
            rtol (float): relative tolerance.

        Raises:
            AssertionError: If the backend under test does not agree with the C++ reference within rtol.
        """
        expected_didv = self._didv_cpp(tip_positions=tip_positions, tip_orb=tip_orb)
        actual_didv = self._didv_generic(tip_positions=tip_positions, tip_orb=tip_orb, backend=didv_backend)

        np.testing.assert_allclose(actual_didv, expected_didv, rtol=rtol)

    def _didv_cpp(self, tip_positions, tip_orb):
        return self._didv_generic(tip_positions=tip_positions, tip_orb=tip_orb, backend=ProbeSTM.dIdV_sp_sp)

    @pytest.fixture(scope="class")
    @abstractmethod
    def all_tip_positions_batch(self):
        pass

    def pytest_generate_tests(self, metafunc):
        if "tip_position" in metafunc.fixturenames:
            tip_positions = list(self._sample_tip_positions(self._NR_SAMPLED_TIP_POSITIONS)
                                 .reshape(self._NR_SAMPLED_TIP_POSITIONS, 1, 1, 1, 3))
            metafunc.parametrize("tip_position", tip_positions,
                                 ids=map(lambda x: str(x.reshape(3)), tip_positions))

        if "tip_orb" in metafunc.fixturenames:
            metafunc.parametrize("tip_orb", self._TIP_ORBITALS, ids=self._TIP_ORBITAL_LABELS)

        if "didv_backend" in metafunc.fixturenames and \
           "rtol" in metafunc.fixturenames:
            metafunc.parametrize(("didv_backend", "rtol"),
                                 (
                                         (ProbeSTMOpenCLParallel.didv,   2e-3),
                                         (ProbeSTMOpenCLSequential.didv, 2e-3),
                                         (ProbeStmNumpy.didv,                           3e-3),
                                         (ProbeStmPytorch.didv,                         3e-3),
                                 ),
                                 ids=[
                                     "OpenCL parallel-rtol=2e-3",
                                     "OpenCL sequential-rtol=2e-3",
                                     "NumPy-rtol=3e-3",
                                     "PyTorch-rtol=3e-3",
                                 ])

    @abstractmethod
    def _sample_tip_positions(self, n_samples: int) -> np.ndarray:
        """Sample n_samples tip positions.

        Args:
            n_samples (int): number of tip positions to sample

        Returns:
            Array with sampled tip positions, shape: (n_samples, 3)
        """
        pass


class _TestDidvBackendConsistencySp(_TestDidvBackendConsistency, ABC):
    """Verify that differential conductance (dI/dV) for sp sample orbitals
        from the C++ backend matches alternative backends."""
    _TIP_ORBITALS = [
        _TestDidvBackendConsistency._S_TIP_ORB,
        _TestDidvBackendConsistency._PZ_TIP_ORB,
        _TestDidvBackendConsistency._PXY_TIP_ORB,
        _TestDidvBackendConsistency._CO_TIP_ORB,
        _TestDidvBackendConsistency._DZ2_TIP_ORB,
        _TestDidvBackendConsistency._DXZYZ_TIP_ORB,
    ]
    _TIP_ORBITAL_LABELS = [
        "s",
        "pz",
        "pxy",
        "CO",
        "dz2",
        "dxzyz",
    ]
    _SAMPLE_ORBITAL_COUNT = 4
    _EIGENENERGIES = np.asarray([-.04, -.03, -.02, -.01])
    _LCAO_KS_COEFFICIENTS_DIAG = [1., 1., 1., 1.]
    _LCAO_KS_COEFFICIENTS = np.diag(_LCAO_KS_COEFFICIENTS_DIAG)


class _TestDidvBackendConsistencySpd(_TestDidvBackendConsistency, ABC):
    """Verify that differential conductance (dI/dV) for spd sample orbitals
    from the C++ backend matches alternative backends."""

    _TIP_ORBITALS = [
        _TestDidvBackendConsistency._S_TIP_ORB,
        _TestDidvBackendConsistency._PZ_TIP_ORB,
        _TestDidvBackendConsistency._PXY_TIP_ORB,
        _TestDidvBackendConsistency._CO_TIP_ORB,
    ]
    _TIP_ORBITAL_LABELS = [
        "s",
        "pz",
        "pxy",
        "CO",
    ]
    _SAMPLE_ORBITAL_COUNT = 9
    _EIGENENERGIES = np.asarray([-.04, -.03, -.02, -.01, 0., .01, .02, .03, .04])
    _LCAO_KS_COEFFICIENTS_DIAG = [1., 1., 1., 1., 0.2, 0.2, 0.2, 0.2, 0.2]
    _LCAO_KS_COEFFICIENTS = np.diag(_LCAO_KS_COEFFICIENTS_DIAG)
