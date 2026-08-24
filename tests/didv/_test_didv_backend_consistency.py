from abc import ABC, abstractmethod
from functools import partial
from typing import Callable, List

import numpy as np
import pytest

from pyPPSTM import ProbeSTM
from pyPPSTM.probe_stm_numpy import ProbeStmNumpy
from pyPPSTM.probe_stm_opencl import ProbeSTMOpenCLParallel, ProbeSTMOpenCLSequential
from pyPPSTM.probe_stm_pytorch import ProbeStmPytorch
from tests.didv import _test_didv


class _TestDidvBackendConsistency(_test_didv._TestDidv, ABC):
    """Shared test scaffold for dI/dV backend consistency checks.

    Subclasses provide concrete fixtures for tip positions and sample orbitals.
    This base class parametrizes backend comparisons against the C++ reference
    implementation."""
    _TIP_ORBITALS: List[np.ndarray]
    _TIP_ORBITAL_LABELS: List[str]

    _NR_SAMPLED_TIP_POSITIONS = 10

    _S_TIP_ORB     = np.asarray([ 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    _PY_TIP_ORB    = np.asarray([ 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    _PZ_TIP_ORB    = np.asarray([ 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    _PX_TIP_ORB    = np.asarray([ 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    _PXY_TIP_ORB   = np.asarray([ 0.0, 0.5, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0])
    _CO_TIP_ORB    = np.asarray([0.15, 0.5, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0])
    _DYZ_TIP_ORB   = np.asarray([ 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0])
    _DZ2_TIP_ORB   = np.asarray([ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0])
    _DXZ_TIP_ORB   = np.asarray([ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0])
    _DXZYZ_TIP_ORB = np.asarray([ 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.5, 0.0])

    _BACKENDS_NAME_FN_RTOL = (
        ("OpenCL parallel",            ProbeSTMOpenCLParallel.didv,                             2e-3),
        ("OpenCL sequential",          ProbeSTMOpenCLSequential.didv,                           2e-3),
        ("NumPy (default chunking)",   ProbeStmNumpy.didv,                                      4e-3),
        ("PyTorch (default chunking)", ProbeStmPytorch.didv,                                    3e-3),
        ("NumPy (no chunking)",        partial(ProbeStmNumpy.didv,    n_tip_position_chunks=1), 4e-3),
        ("PyTorch (no chunking)",      partial(ProbeStmPytorch.didv,  n_tip_position_chunks=1), 3e-3),
        ("NumPy (2 chunks)",           partial(ProbeStmNumpy.didv,    n_tip_position_chunks=2), 4e-3),
        ("PyTorch (2 chunks)",         partial(ProbeStmPytorch.didv,  n_tip_position_chunks=2), 3e-3),
    )

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
                                                                     tip_positions_batch: np.ndarray,
                                                                     tip_orb: np.ndarray,
                                                                     didv_backend: Callable,
                                                                     rtol: float):
        """Verify that dI/dV from the C++ backend matches other backends for a batch of tip positions."""
        self._test_didv_cpp_matches_other_backends(tip_positions_batch,
                                                   tip_orb,
                                                   didv_backend,
                                                   rtol)

    def _test_didv_cpp_matches_other_backends(self,
                                              tip_positions: np.ndarray,
                                              tip_orb: np.ndarray,
                                              didv_backend: Callable,
                                              rtol: float):
        """Compare a backend's dI/dV output against the C++ reference.

        Args:
            tip_positions: Tip positions, shape: (n_z, n_y, n_x, 3).
            tip_orb: Tip-orbital coefficients, shape: (9,).
            didv_backend: Backend implementation under test.
            rtol: Relative tolerance.
        """
        expected_didv = self._didv_cpp(tip_positions=tip_positions, tip_orb=tip_orb)
        actual_didv = self._didv_generic(tip_positions=tip_positions, tip_orb=tip_orb, backend=didv_backend)

        np.testing.assert_allclose(actual_didv, expected_didv, rtol=rtol)

    def _didv_cpp(self, tip_positions, tip_orb):
        return self._didv_generic(tip_positions=tip_positions, tip_orb=tip_orb, backend=ProbeSTM.dIdV_sp_sp)

    @pytest.fixture(scope="class")
    @abstractmethod
    def tip_positions_batch(self):
        """Batch of tip positions."""
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
                                 ((fn, rtol) for _, fn, rtol in self._BACKENDS_NAME_FN_RTOL),
                                 ids=(name for name, _, _ in self._BACKENDS_NAME_FN_RTOL))

    @abstractmethod
    def _sample_tip_positions(self, n_samples: int) -> np.ndarray:
        """Sample tip positions.

        Args:
            n_samples: Number of tip positions to draw.

        Returns:
            Tip positions, shape: (n_samples, 3).
        """
        pass


class _TestDidvBackendConsistencySpSampleOrbitals(_TestDidvBackendConsistency, ABC):
    """Tests dI/dV backend consistency for s and p sample orbitals."""
    _TIP_ORBITALS = [
        _TestDidvBackendConsistency._S_TIP_ORB,
        _TestDidvBackendConsistency._PY_TIP_ORB,
        _TestDidvBackendConsistency._PZ_TIP_ORB,
        _TestDidvBackendConsistency._PX_TIP_ORB,
        _TestDidvBackendConsistency._PXY_TIP_ORB,
        _TestDidvBackendConsistency._CO_TIP_ORB,
        _TestDidvBackendConsistency._DYZ_TIP_ORB,
        _TestDidvBackendConsistency._DZ2_TIP_ORB,
        _TestDidvBackendConsistency._DXZ_TIP_ORB,
        _TestDidvBackendConsistency._DXZYZ_TIP_ORB,
    ]
    _TIP_ORBITAL_LABELS = [
        "s",
        "py",
        "pz",
        "px",
        "pxy",
        "CO",
        "dyz",
        "dz2",
        "dxz",
        "dxzyz",
    ]
    _SAMPLE_ORBITAL_COUNT = 4


class _TestDidvBackendConsistencySpdSampleOrbitals(_TestDidvBackendConsistency, ABC):
    """Tests dI/dV backend consistency for s, p and d sample orbitals."""
    _TIP_ORBITALS = [
        _TestDidvBackendConsistency._S_TIP_ORB,
        _TestDidvBackendConsistency._PY_TIP_ORB,
        _TestDidvBackendConsistency._PZ_TIP_ORB,
        _TestDidvBackendConsistency._PX_TIP_ORB,
        _TestDidvBackendConsistency._PXY_TIP_ORB,
        _TestDidvBackendConsistency._CO_TIP_ORB,
    ]
    _TIP_ORBITAL_LABELS = [
        "s",
        "py",
        "pz",
        "px",
        "pxy",
        "CO",
    ]
    _SAMPLE_ORBITAL_COUNT = 9


class _TestDidvBackendConsistencySpSampleOrbitalsOrbital(_TestDidvBackendConsistencySpSampleOrbitals, ABC):
    """Tests dI/dV backend consistency for s and p sample orbitals.
    Inspired by examples/orbitals_tests."""
    _V = -.04  # voltage

    _EIGENENERGIES = np.asarray([-.04, -.03, -.02, -.01])
    _LCAO_COEFFICIENTS = np.diag([1., 1., 1., 1.])
    _SAMPLE_ATOMS_POSITIONS = np.asarray([[0., 0., 0.]])


class _TestDidvBackendConsistencySpdSampleOrbitalsOrbital(_TestDidvBackendConsistencySpdSampleOrbitals, ABC):
    """Tests dI/dV backend consistency for s, p and d sample orbitals.
    Inspired by examples/orbitals_tests."""
    _V = -.04  # voltage

    _EIGENENERGIES = np.asarray([-.04, -.03, -.02, -.01, 0., .01, .02, .03, .04])
    _LCAO_COEFFICIENTS = np.diag([1., 1., 1., 1., 0.2, 0.2, 0.2, 0.2, 0.2])
    _SAMPLE_ATOMS_POSITIONS = np.asarray([[0., 0., 0.]])
