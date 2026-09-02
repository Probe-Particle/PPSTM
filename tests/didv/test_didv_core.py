from abc import ABC
from functools import partial
from typing import Callable

import numpy as np
import pytest

from ppstm.probe_stm_opencl import ProbeSTMOpenCLParallel, ProbeSTMOpenCLSequential
from ppstm.probe_stm_pytorch import ProbeStmPytorch
from ppstm.probe_stm_numpy import ProbeStmNumpy

# Register assert rewrite BEFORE import to enable better error messages in inherited tests
pytest.register_assert_rewrite("tests.didv._test_didv")
from tests.didv import _test_didv


class _TestDidvCore(_test_didv._TestDidv, ABC):
    """Verify core properties of conductances (dI/dV) computation."""
    _SINGLE_TIP_POSITION: np.ndarray
    _DIFFERENT_TIP_POSITIONS: np.ndarray

    _TIP = np.asarray([1.0] + [0.0]*8) # s tip orbital
    _SAMPLE_ORBITAL_COUNT = 9 # spd = 9 orbitals

    _SAME_TIP_POSITIONS_BATCH_SHAPES = (
        (1, 1, 1),
        (2, 1, 1),
        (1, 3, 1),
        (1, 1, 4),
        (2, 3, 1),
        (1, 3, 4),
        (2, 1, 4),
        (2, 3, 4),
    )

    _ATOL_UB = 1e-7

    _BACKENDS_NAME_FN = (
        ("OpenCL parallel",            ProbeSTMOpenCLParallel.didv),
        ("OpenCL sequential",          ProbeSTMOpenCLSequential.didv),
        ("NumPy (default chunking)",   ProbeStmNumpy.didv),
        ("PyTorch (default chunking)", ProbeStmPytorch.didv),
        ("NumPy (no chunking)",        partial(ProbeStmNumpy.didv,    n_tip_position_chunks=1)),
        ("PyTorch (no chunking)",      partial(ProbeStmPytorch.didv,  n_tip_position_chunks=1)),
        ("NumPy (2 chunks)",           partial(ProbeStmNumpy.didv,    n_tip_position_chunks=2)),
        ("PyTorch (2 chunks)",         partial(ProbeStmPytorch.didv,  n_tip_position_chunks=2)),
    )

    def test_didv_returns_one_conductance_per_tip_position(self,
                                                           single_tip_position: np.ndarray,
                                                           x_len: int, y_len: int, z_len: int,
                                                           didv_backend: Callable):
        """Ensure the dI/dV result contains one conductance value for every tip position
        by matching the spatial array shape."""
        tip_positions = np.tile(single_tip_position, (x_len, y_len, z_len, 1))

        out = self._didv(tip_positions=tip_positions,
                         backend=didv_backend)

        assert out.shape == tip_positions.shape[:3]

    def test_didv_can_differ_for_distinct_tip_positions(self,
                                                        different_tip_positions: np.ndarray,
                                                        didv_backend: Callable):
        """Verify that dI/dV (differential conductance) can differ for distinct tip positions."""
        out = self._didv(tip_positions=different_tip_positions,
                         backend=didv_backend)

        assert out.max() - out.min() > 0.

    def test_didv_deterministic_for_single_tip_position(self,
                                                        single_tip_position: np.ndarray,
                                                        didv_backend: Callable):
        """Verify that dI/dV (differential conductance) is deterministic for a single tip position:
        repeated calls with single tip position produce numerically identical results."""
        out_singles = np.vstack([self._didv(tip_positions=single_tip_position,
                                            backend=didv_backend)
                                 for _ in range(2)])

        assert out_singles.max() - out_singles.min() < self._ATOL_UB

    def test_didv_batch_consistent_for_same_tip_positions(self,
                                                          single_tip_position: np.ndarray,
                                                          x_len: int, y_len: int, z_len: int,
                                                          didv_backend: Callable):
        """Verify that conductances (dI/dV) are constant across the batch
        when all tip positions in the batch are the same."""
        repeated_tip_positions = np.tile(single_tip_position, (x_len, y_len, z_len, 1))
        out = self._didv(tip_positions=repeated_tip_positions,
                         backend=didv_backend)

        assert out.max() - out.min() < self._ATOL_UB

    def test_didv_batch_matches_individual_for_each_tip_position(self,
                                                                 single_tip_position: np.ndarray,
                                                                 different_tip_positions: np.ndarray,
                                                                 are_tip_positions_same: bool,
                                                                 x_len: int, y_len: int, z_len: int,
                                                                 didv_backend: Callable):
        """Verify that computing dI/dV for a batch of different tip positions
        produces the same results as computing dI/dV for each position individually."""
        if are_tip_positions_same:
            tip_positions = np.tile(single_tip_position, (x_len, y_len, z_len, 1))
        else:
            tip_positions = different_tip_positions

        out_multi = self._didv(tip_positions=tip_positions,
                               backend=didv_backend)

        out_singles = np.vstack([self._didv(tip_positions=tip_positions[i:i+1,j:j+1,k:k+1],
                                            backend=didv_backend)
                                 for i in range(x_len) for j in range(y_len) for k in range(z_len)])

        np.testing.assert_equal(out_multi, out_singles.reshape(x_len, y_len, z_len))

    @pytest.fixture(scope="class")
    def single_tip_position(self):
        """A single tip position as a (1,1,1,3) shaped array."""
        yield self._SINGLE_TIP_POSITION

    @pytest.fixture(scope="class")
    def different_tip_positions(self):
        """Multiple tip positions with varying coordinates, shape (2,1,1,3)."""
        yield self._DIFFERENT_TIP_POSITIONS

    def pytest_generate_tests(self, metafunc):
        if "x_len" in metafunc.fixturenames and \
           "y_len" in metafunc.fixturenames and \
           "z_len" in metafunc.fixturenames:
            if "are_tip_positions_same" in metafunc.fixturenames:
                metafunc.parametrize(("are_tip_positions_same", "x_len", "y_len", "z_len"),
                                     (
                                         (False, *self._DIFFERENT_TIP_POSITIONS.shape[:3]),
                                         *[(True, x_len, y_len, z_len)
                                           for x_len, y_len, z_len in self._SAME_TIP_POSITIONS_BATCH_SHAPES],
                                     ),
                                     ids=(
                                         f"vary_tip_pos-("
                                         f"{self._DIFFERENT_TIP_POSITIONS.shape[0]},"
                                         f"{self._DIFFERENT_TIP_POSITIONS.shape[1]},"
                                         f"{self._DIFFERENT_TIP_POSITIONS.shape[2]})",
                                         *[f"same_tip_pos-({x_len},{y_len},{z_len})"
                                           for x_len, y_len, z_len in self._SAME_TIP_POSITIONS_BATCH_SHAPES],
                                     ))
            else:
                metafunc.parametrize(("x_len", "y_len", "z_len"),
                                     self._SAME_TIP_POSITIONS_BATCH_SHAPES)

        if "didv_backend" in metafunc.fixturenames:
            metafunc.parametrize("didv_backend",
                                 (fn for _, fn in self._BACKENDS_NAME_FN),
                                 ids=(name for name, _, in self._BACKENDS_NAME_FN))

    def _didv(self, tip_positions: np.ndarray, backend: Callable) -> np.ndarray:
        """Compute differential conductance using the specified backend.

        Args:
            tip_positions (np.ndarray): array of tip coordinates, shape: (n_z, n_y, n_x, 3)
            backend (Callable): dI/dV implementation

        Returns:
            Array of conductance values, shape: shape (n_z, n_y, n_x)
        """
        return super()._didv_generic(tip_positions=tip_positions,
                                     tip_orb=self._TIP,
                                     backend=backend)


class TestDidvCoreOneSamplePositionDiagonalLCAO(_TestDidvCore):
    """Verify core properties of conductances (dI/dV) computation with a single sample atom position
    and a diagonal LCAO-KS coefficient matrix.
    Inspired by examples/orbitals_tests/orbitals.toml."""
    _V = 0.0     # voltage
    _ETA = 10.0  # eV, Lorentzian width of states in energy scale

    _SINGLE_TIP_POSITION = np.asarray([[[[-5.0, -5.0, -5.0]]]])

    _DIFFERENT_TIP_POSITIONS = np.asarray([
        [[[-5.0, -5.0, -5.0]]],
        [[[-4.8, -4.8, -4.8]]],
    ])

    _EIGENENERGIES = np.asarray([-0.04, -0.03, -0.02, -0.01, 0.0, 0.01, 0.02, 0.03, 0.04])
    _LCAO_COEFFICIENTS = np.diag([1.0, 1.0, 1.0, 1.0, 0.2, 0.2, 0.2, 0.2, 0.2])
    _SAMPLE_ATOMS_POSITIONS = np.asarray([[0.0, 0.0, 0.0]])


class TestDidvCoreTwoSamplePositionsDiagonalLikeLCAO(_TestDidvCore):
    """Verify core properties of conductances (dI/dV) computation with two sample atom positions
    and a LCAO‑KS coefficient matrix formed by horizontally stacking two diagonal matrices.
    Inspired by examples/orbitals_tests/orbitals.toml."""
    _V = 0.0     # voltage
    _ETA = 10.0  # eV, Lorentzian width of states in energy scale

    _SINGLE_TIP_POSITION = np.asarray([[[[-5.0, -5.0, -5.0]]]])

    _DIFFERENT_TIP_POSITIONS = np.asarray([
        [[[-5.0, -5.0, -5.0]]],
        [[[-4.8, -4.8, -4.8]]],
    ])

    _EIGENENERGIES = np.asarray([-0.04, -0.03, -0.02, -0.01, 0.0, 0.01, 0.02, 0.03, 0.04])
    _LCAO_COEFFICIENTS = np.hstack(
        [np.diag([1.0, 1.0, 1.0, 1.0, 0.2, 0.2, 0.2, 0.2, 0.2]),]*2
    )
    _SAMPLE_ATOMS_POSITIONS = np.asarray([[0.0, 0.0, 0.0], ] * 2)


class TestDidvCoreTwoSamplePositionsNonDiagonalLCAO(_TestDidvCore):
    """Verify core properties of conductances (dI/dV) computation with two sample atom positions
    and a non-diagonal LCAO-KS coefficient matrix.
    Inspired by examples/Si_7x7."""
    _V = -2.0   # voltage
    _ETA = 0.1  # eV, Lorentzian width of states in energy scale

    _SINGLE_TIP_POSITION = np.asarray([[[[-2.0, -15.0, 8.7]]]])

    _DIFFERENT_TIP_POSITIONS = np.asarray([
        [[[-2.0, -15.0, 8.7]]],
        [[[-1.8, -15.0, 8.7]]],
    ])

    _EIGENENERGIES = np.asarray([-2.49928599, -2.48934599, -2.48904599, -2.48875599])
    _LCAO_COEFFICIENTS = np.asarray(
        [
            [8.000e-05, -3.251e-02, 1.700e-04, 1.100e-04, -6.900e-04,
             3.760e-04, -8.000e-06, 4.000e-06, 0.000e+00, -1.000e-04,
             1.651e-02, -2.200e-04, 2.829e-02, 3.400e-04, -2.000e-04,
             1.200e-05, -3.200e-04, -5.840e-04],
            [4.400e-04, 6.280e-03, 7.000e-05, 1.790e-03, 5.320e-04,
             -2.300e-04, 1.000e-05, -5.000e-05, -4.400e-05, 4.690e-03,
             2.396e-02, 1.110e-02, 1.703e-02, 2.640e-04, -2.080e-04,
             -5.720e-04, -7.140e-04, -8.060e-04],
            [-5.170e-03, 4.010e-03, -1.326e-02, -1.466e-02, 1.120e-04,
             -1.120e-04, 7.080e-04, -2.380e-04, -2.120e-04, 3.080e-03,
             -8.350e-03, 8.140e-03, -2.214e-02, -6.560e-04, 5.420e-04,
             -4.240e-04, 4.760e-04, 7.180e-04],
            [-1.080e-03, -4.285e-02, -1.630e-03, -3.410e-03, -1.176e-03,
             1.134e-03, 7.200e-05, 2.200e-05, 1.400e-05, -2.610e-03,
             7.740e-03, -5.720e-03, 2.309e-02, 2.700e-04, -3.700e-04,
             2.880e-04, -4.620e-04, -3.240e-04],
        ]
    )
    _SAMPLE_ATOMS_POSITIONS = np.asarray(
        [
            [6.7645490e+00, 1.2580000e-03, 4.6667550e+00],
            [2.0283753e+01, -7.8067410e+00, 4.6567320e+00],
        ]
    )
