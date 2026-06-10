from abc import ABC, abstractmethod
from pathlib import Path
from typing import Optional

import numpy as np
import pytest


class _TestSave(ABC):
    _TIP_TYPE = "fixed"
    _TIP_ORB = "s"
    _ETA = .1
    _WORK_FUNCTION = 5.
    _WORK_FUNCTION_DECAY = .5
    _DIDV = _CURRENT = np.random.randn(1, 1, 1, 1)
    _V = -.5
    _V_MAX = .5
    _DV = .1
    _VOLTAGES = np.arange(_V, _V_MAX + 0.001, _DV)
    _NAMES = [f"{V:.1f}" for V in _VOLTAGES]
    _SCAN_WINDOW = np.asarray([[-5, -5, -5], [5, 5, 5]])
    _LVEC = np.asarray([
        [_SCAN_WINDOW[0, 0], _SCAN_WINDOW[0, 1], _SCAN_WINDOW[0, 2]],
        [_SCAN_WINDOW[1, 0] - _SCAN_WINDOW[0, 0], 0., 0.],
        [0., _SCAN_WINDOW[1, 1] - _SCAN_WINDOW[0, 1], 0.],
        [0., 0., _SCAN_WINDOW[1, 2] - _SCAN_WINDOW[0, 2]]
    ])

    @pytest.fixture(scope="class")
    def current(self, request):
        if request.param:
            return self._CURRENT
        else:
            return None

    @pytest.fixture(scope="class")
    def didv(self, request):
        if request.param:
            return self._DIDV
        else:
            return None

    @pytest.mark.parametrize(("scan_type",     "current", "didv", "expected_files_count"),
                             (("didv",         False,      True,  1),
                              ("dIdV",         False,      True,  1),
                              ("didv-single",  False,      True,  1),
                              ("states",       False,      True,  1),
                              ("STATES",       False,      True,  1),
                              ("STM",          True,       False, 1),
                              ("STM-single",   True,       False, 1),
                              ("v-scan",       True,       True,  2),
                              ("V-scan",       True,       True,  2),
                              ("Voltage-scan", True,       True,  2)),
                             indirect=["current", "didv"])
    def test_all_output_files_count(self, tmp_path: Path, scan_type, current, didv, expected_files_count: int):
        """Verify the number of created output files."""

        self._save(scan_type, current=current, didv=didv, save_dir=tmp_path)

        assert self._count_created_files(tmp_path) == expected_files_count

    @pytest.mark.parametrize(("scan_type",     "current", "didv", "expected_files_count"),
                             (("didv",         False,      True,  1),
                              ("dIdV",         False,      True,  1),
                              ("didv-single",  False,      True,  1),
                              ("states",       False,      True,  1),
                              ("STATES",       False,      True,  1),
                              ("STM",          True,       False, 0),
                              ("STM-single",   True,       False, 0),
                              ("v-scan",       True,       True,  1),
                              ("V-scan",       True,       True,  1),
                              ("Voltage-scan", True,       True,  1)),
                             indirect=["current", "didv"])
    def test_didv_files_count(self, tmp_path: Path, scan_type, current, didv, expected_files_count: int):
        """Verify the number of created output files for dIdV simulations."""

        self._save(scan_type, current=current, didv=didv, save_dir=tmp_path)

        assert self._count_created_files(tmp_path, "didv_*") == expected_files_count

    @pytest.mark.parametrize(("scan_type",     "current", "didv", "expected_files_count"),
                             (("didv",         False,      True,  0),
                              ("dIdV",         False,      True,  0),
                              ("didv-single",  False,      True,  0),
                              ("states",       False,      True,  0),
                              ("STATES",       False,      True,  0),
                              ("STM",          True,       False, 1),
                              ("STM-single",   True,       False, 1),
                              ("v-scan",       True,       True,  1),
                              ("V-scan",       True,       True,  1),
                              ("Voltage-scan", True,       True,  1)),
                             indirect=["current", "didv"])
    def test_stm_files_count(self, tmp_path: Path, scan_type, current, didv, expected_files_count: int):
        """Verify the number of created output files for STM simulations."""

        self._save(scan_type, current=current, didv=didv, save_dir=tmp_path)

        assert self._count_created_files(tmp_path, "STM_*") == expected_files_count

    @pytest.mark.parametrize(("scan_type",     "current", "didv", "expected_work_function_decay"),
                             (("didv",         False,      True,  0.),
                              ("dIdV",         False,      True,  0.),
                              ("didv-single",  False,      True,  0.),
                              ("states",       False,      True,  0.),
                              ("STATES",       False,      True,  0.),
                              ("v-scan",       True,       True,  _WORK_FUNCTION_DECAY),
                              ("V-scan",       True,       True,  _WORK_FUNCTION_DECAY),
                              ("Voltage-scan", True,       True,  _WORK_FUNCTION_DECAY)),
                             indirect=["current", "didv"])
    def test_didv_scan_filename_matches_convention(self, tmp_path: Path, scan_type: str, current, didv, expected_work_function_decay: float):
        """Verify the output filenames for dIdV simulations follow the required format."""

        expected_file_name = self._build_expected_didv_file_name(expected_work_function_decay)

        self._save(scan_type, current=current, didv=didv, save_dir=tmp_path)

        actual_file_name = self._find_one_created_file_name(tmp_path, pattern="didv_*")

        assert actual_file_name == expected_file_name

    @pytest.mark.parametrize(("scan_type",     "current", "didv", "expected_work_function_decay"),
                             (("STM",          True,      False,   _WORK_FUNCTION_DECAY),
                              ("STM-single",   True,      False,   _WORK_FUNCTION_DECAY),
                              ("v-scan",       True,      True,    _WORK_FUNCTION_DECAY),
                              ("V-scan",       True,      True,    _WORK_FUNCTION_DECAY),
                              ("Voltage-scan", True,      True,    _WORK_FUNCTION_DECAY)),
                             indirect=["current", "didv"])
    def test_stm_scan_filename_matches_convention(self, tmp_path: Path, scan_type: str, current, didv, expected_work_function_decay: float):
        """Verify the output filenames for STM simulations follow the required format."""

        expected_file_name = self._build_expected_stm_file_name(expected_work_function_decay)

        self._save(scan_type, current=current, didv=didv, save_dir=tmp_path)

        actual_file_name = self._find_one_created_file_name(tmp_path, pattern="STM_*")

        assert actual_file_name == expected_file_name

    @abstractmethod
    def _save(self, scan_type: str, save_dir: Path, current: Optional = None, didv: Optional = None):
        pass

    @abstractmethod
    def _build_expected_didv_file_name(self, work_function_decay: float) -> str:
        pass

    @abstractmethod
    def _build_expected_stm_file_name(self, expected_work_function_decay: float) -> str:
        pass

    @classmethod
    def _config(cls, scan_type: str):
        return {
            "scan": {
                "tip_type": cls._TIP_TYPE,
                "tip_orb": cls._TIP_ORB,
                "eta": cls._ETA,
                "scan_type": scan_type,
            },
            "advanced": {
                "work_function": cls._WORK_FUNCTION,
                "work_function_decay": cls._WORK_FUNCTION_DECAY,
            }
        }

    @staticmethod
    def _count_created_files(directory: Path, pattern: str = '*') -> int:
        return len(list(directory.glob(pattern)))

    @staticmethod
    def _find_one_created_file_name(directory: Path, pattern: str = '*') -> str:
        return next(iter(directory.glob(pattern))).name
