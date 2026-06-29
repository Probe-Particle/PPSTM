from abc import ABC

import pytest

pytest.register_assert_rewrite("tests.io._test_save")
from tests.io._test_save import _TestSave


class _TestSavePlot(_TestSave, ABC):
    def _build_expected_didv_file_name(self, work_function_decay: float) -> str:
        return f'didv_{self._NAMES[0]}_tip_{self._TIP_TYPE}-{self._TIP_ORB}' \
               f'_WF_{self._WORK_FUNCTION - self._VOLTAGES[0] * work_function_decay}' \
               f'_eta_{self._ETA:.7f}_{0:03d}.{self._FILENAME_EXTENSION}'

    def _build_expected_stm_file_name(self, expected_work_function_decay: float) -> str:
        return f'STM_{self._NAMES[0]}_tip_{self._TIP_TYPE}-{self._TIP_ORB}' \
               f'_WF_{self._WORK_FUNCTION:.1f}_WF_decay_{expected_work_function_decay:.1f}' \
               f'_eta_{self._ETA:.7f}_{0:03d}.{self._FILENAME_EXTENSION}'
