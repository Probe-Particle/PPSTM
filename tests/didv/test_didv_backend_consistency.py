from abc import ABC

import numpy as np
import pytest

# Register assert rewrite BEFORE import to enable better error messages in inherited tests
pytest.register_assert_rewrite("tests.didv._test_didv_backend_consistency")
from tests.didv import _test_didv_backend_consistency


class _TestDidvBackendConsistency(_test_didv_backend_consistency._TestDidvBackendConsistency, ABC):
    """Unit tests for dI/dV backend consistency with synthetic tip positions.

    This test class generates a dense grid of tip positions across a 10×10×10 Ångström region
    (51 positions per axis with 0.2 Ångström spacing) and verifies that the C++ backend
    produces results consistent with OpenCL and NumPy implementations.
    """
    _LOWEST_TIP_POSITION  = -5.0  # Symmetric around origin
    _HIGHEST_TIP_POSITION =  5.0  # Symmetric around origin
    _TIP_POSITION_STEP    =  0.2  # Results in 51×51×51 = 132,651 positions

    _NR_TIP_POSITIONS_SINGLE_AXIS = int((_HIGHEST_TIP_POSITION - _LOWEST_TIP_POSITION) / _TIP_POSITION_STEP) + 1  # 51 positions per axis: (5-(-5))/0.2 + 1 = 51

    # 1D array of tip positions from -5 to 5 Å with 0.2 Å spacing.
    _TIP_POSITIONS_SINGLE_AXIS = np.linspace(start=_LOWEST_TIP_POSITION,
                                             stop=_HIGHEST_TIP_POSITION,
                                             num=_NR_TIP_POSITIONS_SINGLE_AXIS)

    @pytest.fixture(scope="class")
    def all_tip_positions_batch(self):
        """All synthetic tip positions as a grid."""
        yield self._generate_tip_position_grid()

    def _generate_tip_position_grid(self) -> np.ndarray:
        """Generate a 3D grid of tip positions as a Cartesian product of 3 repeated 1D arrays of tip positions."""
        tip_positions_x, tip_positions_y, tip_positions_z = np.meshgrid(
            *((self._TIP_POSITIONS_SINGLE_AXIS,) * 3),
            indexing="ij",  # ensures Cartesian (x,y,z) ordering
        )
        return np.stack((tip_positions_x, tip_positions_y, tip_positions_z), axis=-1)

    def _sample_tip_positions(self, n_samples: int) -> np.ndarray:
        """Sample n_samples synthetic tip positions.

        Args:
            n_samples (int): number of tip positions to sample

        Returns:
            Array with sampled tip positions, shape: (n_samples, 3)

        Raises:
            ValueError: If n_samples exceeds the total number of available positions
        """
        max_positions = self._NR_TIP_POSITIONS_SINGLE_AXIS ** 3

        if n_samples > max_positions:
            raise ValueError(f"Cannot sample {n_samples} positions from {max_positions} available")

        tip_positions_indices = self._make_rng().choice(range(max_positions), size=n_samples, replace=False)
        x_idx, y_idx, z_idx = np.unravel_index(
            tip_positions_indices,
            shape=(self._NR_TIP_POSITIONS_SINGLE_AXIS,) * 3
        )
        return np.column_stack((
            self._TIP_POSITIONS_SINGLE_AXIS[x_idx],
            self._TIP_POSITIONS_SINGLE_AXIS[y_idx],
            self._TIP_POSITIONS_SINGLE_AXIS[z_idx],
        ))


class TestDidvBackendConsistencySp(
    _TestDidvBackendConsistency,
    _test_didv_backend_consistency._TestDidvBackendConsistencySp):
    """Tests dI/dV backend consistency with synthetic tip positions
    for sp sample orbitals (s and p, total 4 eigenstates)."""
    pass


class TestDidvBackendConsistencySpd(
    _TestDidvBackendConsistency,
    _test_didv_backend_consistency._TestDidvBackendConsistencySpd):
    """Tests dI/dV backend consistency with synthetic tip positions
    for spd sample orbitals (s, p, and d, total 9 eigenstates)."""
    pass
