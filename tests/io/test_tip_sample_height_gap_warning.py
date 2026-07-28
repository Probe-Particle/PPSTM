import logging
from unittest.mock import patch, MagicMock

import numpy as np

import ppstm_run


class TestTipSampleHeightGapWarning:
    _RNG = np.random.default_rng(42)

    _SAMPLE_ATOM_POSITIONS = _RNG.integers(low=1, high=5, size=(2, 3))

    _MIN_TIP_SAMPLE_HEIGHT_GAP = 2.0  # Angstroms
    _SAFE_TIP_HEIGHT   = np.max(_SAMPLE_ATOM_POSITIONS[...,2]) + _MIN_TIP_SAMPLE_HEIGHT_GAP
    _UNSAFE_TIP_HEIGHT = np.max(_SAMPLE_ATOM_POSITIONS[...,2]) + _MIN_TIP_SAMPLE_HEIGHT_GAP - 0.1

    _SAFE_TIP_POSITIONS = np.hstack((
        _RNG.integers(low=1, high=5, size=(8, 2)),
        _SAFE_TIP_HEIGHT + np.arange(0, 8).reshape(8, 1),
    )).reshape(2, 2, 2, 3)

    _UNSAFE_TIP_POSITIONS = _SAFE_TIP_POSITIONS.copy()
    _UNSAFE_TIP_POSITIONS[0, 0, 0, 2] = _UNSAFE_TIP_HEIGHT

    _WARNING_MESSAGE = "Detected minimum tip height below the maximum sample height " \
                       f"within {_MIN_TIP_SAMPLE_HEIGHT_GAP} Å"

    def test_tip_below_sample_height_raise_warning(self, caplog):
        self._run_ppstm(self._UNSAFE_TIP_POSITIONS)
        assert (ppstm_run.logger.name, logging.WARNING, self._WARNING_MESSAGE) in caplog.record_tuples

    def test_tip_above_sample_height_no_warning(self, caplog):
        self._run_ppstm(self._SAFE_TIP_POSITIONS)
        assert (ppstm_run.logger.name, logging.WARNING, self._WARNING_MESSAGE) not in caplog.record_tuples

    def _run_ppstm(self, tip_positions):
        with patch("pyPPSTM.STMutils.get_tip_positions", return_value=(tip_positions, None, None, None, None)), \
             patch("pyPPSTM.ReadSTM.read_dft", return_value=(None, None, self._SAMPLE_ATOM_POSITIONS)), \
             patch("pyPPSTM.STMutils.get_tip_coefficients"), \
             patch("pyPPSTM.STMutils.run_stm_scan", return_value=(None, None)), \
             patch("pyPPSTM.visualization.get_voltages_and_names", return_value=(None, None)), \
             patch("pyPPSTM.visualization.plot_png"), \
             patch("pyPPSTM.visualization.plot_wsxm"), \
             patch("pyPPSTM.visualization.save_xsf"), \
             patch("pyPPSTM.visualization.save_npz"):
            ppstm_run.run_simulation(config=MagicMock())
