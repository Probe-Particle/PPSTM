from abc import ABC

import numpy as np
import pytest

from pyPPSTM import STMutils

# Register assert rewrite BEFORE import to enable better error messages in inherited tests
pytest.register_assert_rewrite("tests.didv._test_didv_backend_consistency")
from tests.didv import _test_didv_backend_consistency


class _TestDidvBackendConsistencyWithRealTipPositions(_test_didv_backend_consistency._TestDidvBackendConsistency, ABC):
    """Integration tests for dI/dV backend consistency with real tip positions.

    This test class delegates to the official STMutils.get_tip_positions()
    the generation of a dense grid of tip positions across a 10×10×10 Ångström region
    (51 positions per axis with 0.2 Ångström spacing) and verifies that the C++ backend
    produces results consistent with OpenCL and NumPy implementations.
    """
    _TIP_TYPE = "fixed"

    # Scan parameters: 10 Å × 10 Å × 10 Å volume with 50×50×50 point grid
    _SCAN_WINDOW = [
        [-5, -5, -5],
        [ 5,  5,  5]
    ]  # Angstroms
    _SCAN_DIM = [50, 50, 50]  # Points per dimension

    _PLOT_ATOMS = False

    @pytest.fixture(scope="class")
    def all_tip_positions_batch(self):
        """All real (STMutils-sourced) tip positions."""
        yield self._batch_all_tip_positions()

    def _sample_tip_positions(self, n_samples: int) -> np.ndarray:
        """Sample n_samples tip positions from the real (STMutils) batch.

        Args:
            n_samples (int): number of tip positions to sample

        Returns:
            Array of shape (n_samples, 3) with sampled tip positions

        Raises:
            ValueError: If n_samples exceeds the total number of available positions
        """
        all_tip_positions = self._batch_all_tip_positions().reshape(-1, 3)

        if n_samples > (max_positions := len(all_tip_positions)):
            raise ValueError(f"Cannot sample {n_samples} positions from {max_positions} available")

        tip_positions_indices = self._make_rng().choice(range(max_positions), size=n_samples, replace=False)
        return all_tip_positions[tip_positions_indices]

    def _batch_all_tip_positions(self) -> np.ndarray:
        """All real (STMutils-sourced) tip positions.

        Returns:
            Array with all tip positions, shape: (n_z, n_y, n_x, 3)
        """
        tip_r, tip_r0, lvec, extent, atomic_head_or_info = STMutils.get_tip_positions({
            "scan": {
                "tip_type": self._TIP_TYPE,
                "scan_window": self._SCAN_WINDOW,
                "scan_dim": self._SCAN_DIM,
            },
            "output": {
                "plot_atoms": self._PLOT_ATOMS
            },
        })
        return tip_r


class TestDidvBackendConsistencyWithRealTipPositionsSp(
    _TestDidvBackendConsistencyWithRealTipPositions,
    _test_didv_backend_consistency._TestDidvBackendConsistencySp):
    """Tests dI/dV backend consistency with real (STMutils-sourced) tip positions
    for sp sample orbitals (s and p, total 4 eigenstates)."""
    pass


class TestDidvBackendConsistencyWithRealTipPositionsSpd(
    _TestDidvBackendConsistencyWithRealTipPositions,
    _test_didv_backend_consistency._TestDidvBackendConsistencySpd):
    """Tests dI/dV backend consistency with real (STMutils-sourced) tip positions
    for spd sample orbitals (s, p, and d, total 9 eigenstates)."""
    pass
