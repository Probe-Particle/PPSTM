from pathlib import Path
from typing import Optional

import numpy as np
import pytest

from pyPPSTM import visualization

pytest.register_assert_rewrite("tests.io._test_plot")
from tests.io._test_save_plot import _TestSavePlot


class TestSavePlotWsxm(_TestSavePlot):
    _DIDV = _CURRENT = np.random.randn(1, 1, 4, 1)
    _SCAN_DIM = [1, 1, 1]
    _DX = (_TestSavePlot._SCAN_WINDOW[1, 0] - _TestSavePlot._SCAN_WINDOW[0, 0]) / _SCAN_DIM[0]
    _DY = (_TestSavePlot._SCAN_WINDOW[1, 1] - _TestSavePlot._SCAN_WINDOW[0, 1]) / _SCAN_DIM[1]
    _DZ = (_TestSavePlot._SCAN_WINDOW[1, 2] - _TestSavePlot._SCAN_WINDOW[0, 2]) / _SCAN_DIM[2]
    _TIP_R0 = np.mgrid[
        _TestSavePlot._SCAN_WINDOW[0, 0]:_TestSavePlot._SCAN_WINDOW[1, 0] + 0.0001:_DX,
        _TestSavePlot._SCAN_WINDOW[0, 1]:_TestSavePlot._SCAN_WINDOW[1, 1] + 0.0001:_DY,
        _TestSavePlot._SCAN_WINDOW[0, 2]:_TestSavePlot._SCAN_WINDOW[1, 2] + 0.0001:_DZ
    ].transpose()
    _FILENAME_EXTENSION = "wsxm"

    def _save(self, scan_type: str, save_dir: Path, current: Optional = None, didv: Optional = None):
        visualization.plot_wsxm(self._config(scan_type),
                                current=current,
                                didv=didv,
                                voltages=self._VOLTAGES,
                                names=self._NAMES,
                                tip_r0=self._TIP_R0,
                                save_dir=save_dir)
