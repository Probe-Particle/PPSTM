from pathlib import Path
from typing import Optional

import pytest

from pyPPSTM import visualization

pytest.register_assert_rewrite("tests.io._test_plot")
from tests.io._test_save_plot import _TestSavePlot


class TestSavePlotPng(_TestSavePlot):
    _FILENAME_EXTENSION = "png"

    def _save(self, scan_type: str, save_dir: Path, current: Optional = None, didv: Optional = None):
        visualization.plot_png(self._config(scan_type),
                               current=current,
                               didv=didv,
                               voltages=self._VOLTAGES,
                               names=self._NAMES,
                               lvec=self._LVEC,
                               extent=None,
                               geom_plot=None,
                               save_dir=save_dir)
