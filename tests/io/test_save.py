from abc import ABC, abstractmethod
from pathlib import Path
from typing import Optional, List, Tuple

import numpy as np
import pytest

from pyPPSTM import visualization


class _TestSave(ABC):
    _FILENAME_EXTENSION: str
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

    @pytest.fixture(scope="function")
    def current(self, request):
        if request.param:
            yield self._CURRENT
        else:
            yield None

    @pytest.fixture(scope="function")
    def didv(self, request):
        if request.param:
            yield self._DIDV
        else:
            yield None

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
    def test_all_output_files_count(self,
                                    tmp_path: Path,
                                    scan_type,
                                    current,
                                    didv,
                                    expected_files_count: int,
                                    geom_plot: Tuple[List[str], List[float], List[float], List[float], List[float]]|None):
        """Verify the number of created output files."""

        self._save(scan_type, current=current, didv=didv, save_dir=tmp_path, geom_plot=geom_plot)

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
    def test_didv_files_count(self,
                              tmp_path: Path,
                              scan_type,
                              current,
                              didv,
                              expected_files_count: int,
                              geom_plot: Tuple[List[str], List[float], List[float], List[float], List[float]]|None):
        """Verify the number of created output files for dIdV simulations."""

        self._save(scan_type, current=current, didv=didv, save_dir=tmp_path, geom_plot=geom_plot)

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
    def test_stm_files_count(self,
                             tmp_path: Path,
                             scan_type,
                             current,
                             didv,
                             expected_files_count: int,
                             geom_plot: Tuple[List[str], List[float], List[float], List[float], List[float]]|None):
        """Verify the number of created output files for STM simulations."""

        self._save(scan_type, current=current, didv=didv, save_dir=tmp_path, geom_plot=geom_plot)

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
    def test_didv_scan_filename_matches_convention(self,
                                                   tmp_path: Path,
                                                   scan_type: str,
                                                   current,
                                                   didv,
                                                   expected_work_function_decay: float,
                                                   geom_plot: Tuple[List[str], List[float], List[float], List[float], List[float]]|None):
        """Verify the output filenames for dIdV simulations follow the required format."""

        expected_file_name = self._build_expected_didv_file_name(expected_work_function_decay, geom_plot=geom_plot)

        self._save(scan_type, current=current, didv=didv, save_dir=tmp_path, geom_plot=geom_plot)

        actual_file_name = self._find_one_created_file_name(tmp_path, pattern="didv_*")

        assert actual_file_name == expected_file_name

    @pytest.mark.parametrize(("scan_type",     "current", "didv", "expected_work_function_decay"),
                             (("STM",          True,      False,   _WORK_FUNCTION_DECAY),
                              ("STM-single",   True,      False,   _WORK_FUNCTION_DECAY),
                              ("v-scan",       True,      True,    _WORK_FUNCTION_DECAY),
                              ("V-scan",       True,      True,    _WORK_FUNCTION_DECAY),
                              ("Voltage-scan", True,      True,    _WORK_FUNCTION_DECAY)),
                             indirect=["current", "didv"])
    def test_stm_scan_filename_matches_convention(self,
                                                  tmp_path: Path,
                                                  scan_type: str,
                                                  current,
                                                  didv,
                                                  expected_work_function_decay: float,
                                                  geom_plot: Tuple[List[str], List[float], List[float], List[float], List[float]]|None):
        """Verify the output filenames for STM simulations follow the required format."""

        expected_file_name = self._build_expected_stm_file_name(expected_work_function_decay, geom_plot=geom_plot)

        self._save(scan_type, current=current, didv=didv, save_dir=tmp_path, geom_plot=geom_plot)

        actual_file_name = self._find_one_created_file_name(tmp_path, pattern="STM_*")

        assert actual_file_name == expected_file_name

    @abstractmethod
    def _save(self,
              scan_type: str,
              save_dir: Path,
              current: Optional = None,
              didv: Optional = None,
              geom_plot: Tuple[List[str], List[float], List[float], List[float], List[float]]|None = None):
        pass

    @abstractmethod
    def _build_expected_didv_file_name(self,
                                       work_function_decay: float,
                                       geom_plot: Tuple[List[str], List[float], List[float], List[float], List[float]]|None = None) -> str:
        pass

    @abstractmethod
    def _build_expected_stm_file_name(self,
                                      work_function_decay: float,
                                      geom_plot: Tuple[List[str], List[float], List[float], List[float], List[float]]|None = None) -> str:
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

    def pytest_generate_tests(self, metafunc):
        if "geom_plot" in metafunc.fixturenames:
            metafunc.parametrize("geom_plot", [
                None,
                (['C','H'], [11.501016, 3.798985], [6.729982, 11.15007], [5.0, 5.0], [0.0, 0.0]),
            ], ids=(
                "",
                "geom",
            ))

class _TestSaveData(_TestSave, ABC):
    def _build_expected_didv_file_name(self,
                                       work_function_decay: float,
                                       geom_plot: Tuple[List[str], List[float], List[float], List[float], List[float]]|None = None) -> str:
        return f'didv_{self._NAMES[0]}_tip_{self._TIP_TYPE}-{self._TIP_ORB}' \
               f'_WF_{self._WORK_FUNCTION - self._VOLTAGES[0] * work_function_decay}' \
               f'_eta_{self._ETA:.7f}.{self._FILENAME_EXTENSION}'

    def _build_expected_stm_file_name(self,
                                      work_function_decay: float,
                                      geom_plot: Tuple[List[str], List[float], List[float], List[float], List[float]]|None = None) -> str:
        return f'STM_{self._NAMES[0]}_tip_{self._TIP_TYPE}-{self._TIP_ORB}' \
               f'_WF_{self._WORK_FUNCTION:.1f}_WF_decay_{work_function_decay:.1f}' \
               f'_eta_{self._ETA:.7f}.{self._FILENAME_EXTENSION}'


class _TestSavePlot(_TestSave, ABC):
    def _build_expected_didv_file_name(self,
                                       work_function_decay: float,
                                       geom_plot: Tuple[List[str], List[float], List[float], List[float], List[float]]|None = None) -> str:
        return f'didv_{self._NAMES[0]}_tip_{self._TIP_TYPE}-{self._TIP_ORB}' \
               f'_WF_{self._WORK_FUNCTION - self._VOLTAGES[0] * work_function_decay}' \
               f'_eta_{self._ETA:.7f}_{0:03d}.{self._FILENAME_EXTENSION}'

    def _build_expected_stm_file_name(self,
                                      work_function_decay: float,
                                      geom_plot: Tuple[List[str], List[float], List[float], List[float], List[float]]|None = None) -> str:
        return f'STM_{self._NAMES[0]}_tip_{self._TIP_TYPE}-{self._TIP_ORB}' \
               f'_WF_{self._WORK_FUNCTION:.1f}_WF_decay_{work_function_decay:.1f}' \
               f'_eta_{self._ETA:.7f}_{0:03d}.{self._FILENAME_EXTENSION}'


class _TestSaveGeometry(_TestSave, ABC):
    def _build_expected_didv_file_name(self,
                                       work_function_decay: float,
                                       geom_plot: Tuple[List[str], List[float], List[float], List[float], List[float]]|None = None) -> str:
        return self._build_expected_file_name(super()._build_expected_didv_file_name(work_function_decay),
                                              geom_plot=geom_plot)

    def _build_expected_stm_file_name(self,
                                      work_function_decay: float,
                                      geom_plot: Tuple[List[str], List[float], List[float], List[float], List[float]]|None = None) -> str:
        return self._build_expected_file_name(super()._build_expected_stm_file_name(work_function_decay),
                                              geom_plot=geom_plot)

    @classmethod
    def _build_expected_file_name(cls,
                                  file_name: str,
                                  geom_plot: Tuple[List[str], List[float], List[float], List[float], List[float]]|None) -> str:
        return f"{Path(file_name).stem}_geom.{cls._FILENAME_EXTENSION}" if geom_plot is not None else file_name


class TestSaveNpz(_TestSaveData):
    _ATOMIC_INFO_OR_HEAD = (np.zeros((4, 1)), np.zeros((4, 3)))
    _FILENAME_EXTENSION = "npz"

    def _save(self,
              scan_type: str,
              save_dir: Path,
              current: Optional = None,
              didv: Optional = None,
              geom_plot: Tuple[List[str], List[float], List[float], List[float], List[float]]|None = None):
        visualization.save_npz(self._config(scan_type),
                               current=current,
                               didv=didv,
                               voltages=self._VOLTAGES,
                               names=self._NAMES,
                               lvec=self._LVEC,
                               atomic_info_or_head=self._ATOMIC_INFO_OR_HEAD,
                               save_dir=save_dir)


class TestSavePlotPng(_TestSaveGeometry, _TestSavePlot):
    _FILENAME_EXTENSION = "png"

    def _save(self,
              scan_type: str,
              save_dir: Path,
              current: Optional = None,
              didv: Optional = None,
              geom_plot: Tuple[List[str], List[float], List[float], List[float], List[float]]|None = None):
        visualization.plot_png(self._config(scan_type),
                               current=current,
                               didv=didv,
                               voltages=self._VOLTAGES,
                               names=self._NAMES,
                               lvec=self._LVEC,
                               extent=None,
                               geom_plot=geom_plot,
                               save_dir=save_dir)


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

    def _save(self,
              scan_type: str,
              save_dir: Path,
              current: Optional = None,
              didv: Optional = None,
              geom_plot: Tuple[List[str], List[float], List[float], List[float], List[float]]|None = None):
        visualization.plot_wsxm(self._config(scan_type),
                                current=current,
                                didv=didv,
                                voltages=self._VOLTAGES,
                                names=self._NAMES,
                                tip_r0=self._TIP_R0,
                                save_dir=save_dir)


class TestSaveXsf(_TestSaveGeometry, _TestSaveData):
    _FILENAME_EXTENSION = "xsf"

    def _save(self,
              scan_type: str,
              save_dir: Path,
              current: Optional = None,
              didv: Optional = None,
              geom_plot: Tuple[List[str], List[float], List[float], List[float], List[float]]|None = None):
        visualization.save_xsf(self._config(scan_type),
                               current=current,
                               didv=didv,
                               voltages=self._VOLTAGES,
                               names=self._NAMES,
                               lvec=self._LVEC,
                               geom_plot=geom_plot,
                               save_dir=save_dir)
