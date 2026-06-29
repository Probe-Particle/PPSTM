from pathlib import Path
from typing import Optional

import pytest

from pyPPSTM import visualization

pytest.register_assert_rewrite("tests.io._test_save_data")
from tests.io._test_save_data import _TestSaveData


class TestSaveXsf(_TestSaveData):
    _FILENAME_EXTENSION = "xsf"

    def _save(self, scan_type: str, save_dir: Path, current: Optional = None, didv: Optional = None):
        visualization.save_xsf(self._config(scan_type),
                               current=current,
                               didv=didv,
                               voltages=self._VOLTAGES,
                               names=self._NAMES,
                               lvec=self._LVEC,
                               geom_plot=None,
                               save_dir=save_dir)

    def _build_expected_didv_file_name(self, work_function_decay: float) -> str:
        return f'didv_{self._NAMES[0]}_tip_{self._TIP_TYPE}-{self._TIP_ORB}' \
               f'_WF_{self._WORK_FUNCTION - self._VOLTAGES[0] * work_function_decay}' \
               f'_eta_{self._ETA:.7f}.{self._FILENAME_EXTENSION}'

    def _build_expected_stm_file_name(self, expected_work_function_decay: float) -> str:
        return f'STM_{self._NAMES[0]}_tip_{self._TIP_TYPE}-{self._TIP_ORB}' \
               f'_WF_{self._WORK_FUNCTION:.1f}_WF_decay_{expected_work_function_decay:.1f}' \
               f'_eta_{self._ETA:.7f}.{self._FILENAME_EXTENSION}'
