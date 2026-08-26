import logging
from collections import defaultdict
from unittest.mock import patch

import numpy as np

import ppstm_run


@patch('pyPPSTM.STMutils.get_tip_coefficients')
@patch('pyPPSTM.STMutils.run_stm_scan', return_value=(None, None))
@patch('pyPPSTM.visualization.get_voltages_and_names', return_value=(None, None))
class TestTipSampleHeightGapWarning:
    _WARNING_MESSAGE = "Detected minimum tip height below the maximum sample height " \
                       f"within {ppstm_run._MIN_SAMPLE_ATOM_HEIGHTS_GAP} Å"

    @patch('pyPPSTM.ReadSTM.read_dft',
           return_value=(None, None, np.asarray([
               [1, 2, 3],
               [3, 4, 5],
           ])))
    @patch('pyPPSTM.STMutils.get_tip_positions',
           return_value=(
                   np.asarray([
                       [5, 6, 6.999],
                       [7, 8, 7],
                       [9, 1, 8],
                       [1, 3, 8],
                       [1, 4, 9],
                       [1, 6, 9],
                       [7, 1, 10],
                       [3, 8, 10],
                   ]).reshape(2, 2, 2, 3),
                   None,
                   None,
                   None,
                   None,
           ))
    def test_tip_below_sample_height_raise_warning(self, mock1, mock2, mock3, mock4, mock5, caplog):
        self._run_ppstm()
        assert (ppstm_run.logger.name, logging.WARNING, self._WARNING_MESSAGE) in caplog.record_tuples

    @patch('pyPPSTM.ReadSTM.read_dft',
           return_value=(None, None, np.asarray([
               [1, 2, 3],
               [3, 4, 5],
           ])))
    @patch('pyPPSTM.STMutils.get_tip_positions',
           return_value=(
                   np.asarray([
                       [5, 6, 7],
                       [7, 8, 7],
                       [9, 1, 8],
                       [1, 3, 8],
                       [1, 4, 9],
                       [1, 6, 9],
                       [7, 1, 10],
                       [3, 8, 10],
                   ]).reshape(2, 2, 2, 3),
                   None,
                   None,
                   None,
                   None,
           ))
    def test_tip_above_sample_height_no_warning(self, mock1, mock2, mock3, mock4, mock5, caplog):
        self._run_ppstm()
        assert (ppstm_run.logger.name, logging.WARNING, self._WARNING_MESSAGE) not in caplog.record_tuples

    def _run_ppstm(self):
        ppstm_run.main(config=defaultdict(lambda: defaultdict(lambda: None)))
