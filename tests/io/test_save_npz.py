from pathlib import Path
from typing import Optional

import numpy as np
import pytest

from pyPPSTM import visualization

pytest.register_assert_rewrite("tests.io._test_save_data")
from tests.io._test_save_data import _TestSaveData


class TestSaveNpz(_TestSaveData):
    _ATOMIC_INFO_OR_HEAD = (np.zeros((4, 1)), np.zeros((4, 3)))
    _FILENAME_EXTENSION = "npz"

    def _save(self, scan_type: str, save_dir: Path, current: Optional = None, didv: Optional = None):
        visualization.save_npy(self._config(scan_type),
                               current=current,
                               didv=didv,
                               voltages=self._VOLTAGES,
                               names=self._NAMES,
                               lvec=self._LVEC,
                               atomic_info_or_head=self._ATOMIC_INFO_OR_HEAD,
                               save_dir=save_dir)
