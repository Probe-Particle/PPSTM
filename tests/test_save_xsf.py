from pathlib import Path
from typing import Optional

import numpy as np
import pytest

from pyPPSTM import visualization

_TIP_TYPE = "fixed"
_TIP_ORB = "s"
_ETA = .1
_WORK_FUNCTION = 5.
_WORK_FUNCTION_DECAY = .5

_DIDV = _CURRENT = np.random.randn(1, 1, 1, 1)
_V = -.5
_V_MAX = .5
_DV =  .1
_VOLTAGES = np.arange(_V, _V_MAX + 0.001, _DV)
_NAMES = [f"{V:.1f}" for V in _VOLTAGES]
_SCAN_WINDOW = np.asarray([[-5, -5, -5], [5, 5, 5]])
_LVEC = np.asarray([
    [_SCAN_WINDOW[0,0], _SCAN_WINDOW[0,1], _SCAN_WINDOW[0,2]],
    [_SCAN_WINDOW[1,0]-_SCAN_WINDOW[0,0], 0., 0.],
    [0., _SCAN_WINDOW[1,1]-_SCAN_WINDOW[0,1], 0.],
    [0., 0., _SCAN_WINDOW[1,2]-_SCAN_WINDOW[0,2]]
])

_FILENAME_EXTENSION = "xsf"


@pytest.mark.parametrize("scan_type", ("didv", "dIdV", "didv-single", "states", "STATES"))
def test_didv_states_scan_filename_matches_convention(tmp_path: Path, scan_type):
    """Verify XSF filenames for dIdV and states scans follow the required format."""

    expected_file_name = _build_expected_didv_file_name(work_function_decay=0.)

    _save_xsf(scan_type,
              didv=_DIDV,
              save_dir=tmp_path)

    assert len(list(tmp_path.glob("*"))) == 1
    assert next(iter(tmp_path.glob("*"))).name == expected_file_name


@pytest.mark.parametrize("scan_type", ("STM", "STM-single"))
def test_stm_scan_filename_matches_convention(tmp_path: Path, scan_type):
    """Verify XSF filenames for STM scans follow the required format."""

    expected_file_name = _build_expected_stm_file_name()

    _save_xsf(scan_type,
              current=_CURRENT,
              save_dir=tmp_path)

    assert len(list(tmp_path.glob("*"))) == 1
    assert next(iter(tmp_path.glob("*"))).name == expected_file_name


@pytest.mark.parametrize("scan_type", ("v-scan", "V-scan", "Voltage-scan"))
def test_voltage_scan_filename_matches_convention(tmp_path: Path, scan_type):
    """Verify XSF filenames for voltage scans follow the required format."""

    expected_didv_file_name = _build_expected_didv_file_name()
    expected_current_file_name = _build_expected_stm_file_name()

    _save_xsf(scan_type,
              current=_CURRENT,
              didv = _DIDV,
              save_dir=tmp_path)

    assert len(list(tmp_path.glob("*"))) == 2
    assert len(list(tmp_path.glob("didv_*"))) == 1
    assert len(list(tmp_path.glob("STM_*"))) == 1
    assert next(iter(tmp_path.glob("didv_*"))).name == expected_didv_file_name
    assert next(iter(tmp_path.glob("STM_*"))).name == expected_current_file_name


def _build_expected_didv_file_name(work_function_decay: float = _WORK_FUNCTION_DECAY) -> str:
    return f'didv_{_NAMES[0]}_tip_{_TIP_TYPE}-{_TIP_ORB}' \
           f'_WF_{_WORK_FUNCTION - _VOLTAGES[0] * work_function_decay}_eta_{_ETA:.1f}.{_FILENAME_EXTENSION}'


def _build_expected_stm_file_name() -> str:
    return f'STM_{_NAMES[0]}_tip_{_TIP_TYPE}-{_TIP_ORB}' \
           f'_WF_{_WORK_FUNCTION:.1f}_WF_decay_{_WORK_FUNCTION_DECAY:.1f}' \
           f'_eta_{_ETA:.1f}.{_FILENAME_EXTENSION}'

def _save_xsf(scan_type: str, save_dir: Path, current: Optional = None, didv: Optional = None):
    visualization.save_xsf(_config(scan_type),
                           current=current,
                           didv=didv,
                           voltages=_VOLTAGES,
                           names=_NAMES,
                           lvec=_LVEC,
                           geom_plot=None,
                           save_dir=save_dir)


def _config(scan_type: str):
    return {
        "scan": {
            "tip_type": _TIP_TYPE,
            "tip_orb": _TIP_ORB,
            "eta": _ETA,
            "scan_type": scan_type,
        },
        "advanced": {
            "work_function": _WORK_FUNCTION,
            "work_function_decay": _WORK_FUNCTION_DECAY,
        }
    }
