import itertools
from collections import defaultdict
from pathlib import Path
from typing import Iterable, List, Dict, Any

import numpy as np
import toml

from pyPPSTM import DidvBackend
from pyPPSTM.scan_config_builder import ScanConfigBuilder


class TestScanConfigBuilder:
    _TOML_ARGS = {
        "scan_type":     "scan_type",
        "tip_type":      "tip_type",
        "scan_window":   "scan_window",
        "scan_dim":      "scan_dim",
        "v":             "V",
        "v_max":         "V_max",
        "dv":            "dV",
        "eta":           "eta",
        "tip_orb":       "tip_orb",
        "sample_orbs":   "sample_orbs",
        "path":          "path",
        "dft_code":      "dft_code",
        "geometry_file": "geometry_file",
        "pbc":           "pbc",
        "lvs":           "lvs",
        "spin":          "spin",
        "cp2k_name":     "cp2k_name",
        "q":             "Q",
        "k":             "K",
        "data_format":   "data_format",
        "cut_atoms":     "cut_atoms",
        "lower_atoms":   "lower_atoms",
        "lower_coefs":   "lower_coefs",
        "work_function": "work_function",
        "fermi":         "fermi",
        "cut_min":       "cut_min",
        "cut_max":       "cut_max",
        "didv_backend":  "didv_backend",
        'plot_atoms':    "plot_atoms",
        "save_npy":      "NPY",
        "save_png":      "PNG",
        "save_wsxm":     "WSxM",
        "save_xsf":      "XSF",
    }

    _ARGS_GROUP = {
        "scan_type":     "scan",
        "tip_type":      "scan",
        "scan_window":   "scan",
        "scan_dim":      "scan",
        "v":             "scan",
        "v_max":         "scan",
        "dv":            "scan",
        "eta":           "scan",
        "tip_orb":       "scan",
        "sample_orbs":   "scan",
        "path":          "input",
        "dft_code":      "input",
        "geometry_file": "input",
        "pbc":           "input",
        "lvs":           "input",
        "spin":          "input",
        "cp2k_name":     "input",
        "q":             "input",
        "k":             "input",
        "data_format":   "input",
        "cut_atoms":     "advanced",
        "lower_atoms":   "advanced",
        "lower_coefs":   "advanced",
        "work_function": "advanced",
        "fermi":         "advanced",
        "cut_min":       "advanced",
        "cut_max":       "advanced",
        "didv_backend":  "advanced",
        'plot_atoms':    "output",
        "save_npy":      "output",
        "save_png":      "output",
        "save_wsxm":     "output",
        "save_xsf":      "output",
    }

    _TIP_TYPE_CHOICES = [
        "fixed",
        "relaxed",
    ]

    _TIP_ORB_CHOICES = [
        "s",
        "pz",
        "pxy",
        "spxy",
        "5spxy",
        "5spxy",
        "10spxy",
        "CO",
        "dz2",
        "dxzyz"
    ]

    _SAMPLE_ORBS_CHOICES = [
        "sp",
        "spd",
    ]

    _DFT_CODE_CHOICES = [
        "fireball",
        "cp2k",
        "gpaw",
        "aims",
    ]

    _PBC_CHOICES = [
        (0, 0),
        (0.5, 0.5),
        (1, 1),
        (2, 2),
        (3, 3),
    ]

    _DEFAULT_ARG_VALUES = {
        "cut_atoms":   -1,
        "lower_atoms": [],
        "lower_coefs": [],
        "spin":        None,
        "cp2k_name":   "",
        "didv_backend": DidvBackend.CPP,
        "plot_atoms":  False,
        "save_npy":    False,
        "save_png":    False,
        "save_wsxm":   False,
        "save_xsf":    False,
    }

    _BASE_CONFIG_EXAMPLE = {
        "scan_type": "states",
        "v": -0.05,
        "v_max": 0.05,
        "dv": 0.1,
        "eta": 1e-7,
        "path": Path("./"),
        "geometry_file": "input_plot.xyz",
        "work_function": 5,
        "fermi": 0,
        "cut_min": -2.5,
        "cut_max": 2.5,
    }

    _CONFIG_EXAMPLE_BY_TIP_TYPE = {
        "fixed": {
            "tip_type": "fixed",
            "scan_window": np.asarray([[-5, -5, -5], [5, 5, 5]]),
            "scan_dim": [50, 50, 50],
        },
        "relaxed": {
            "tip_type": "relaxed",
            "q": 0,
            "k": 0.24,
            "data_format": "npy",
        },
    }

    _CONFIG_EXAMPLE_BY_DFT_CODE = {
        "fireball": {
            "dft_code": "fireball",
            "lvs": "false",
        },
        "cp2k": {
            "dft_code": "cp2k",
            "cp2k_name": "crazy_mol",
            "lvs": "false",
        },
        "gpaw": {
            "dft_code": "gpaw",
        },
        "aims": {
            "dft_code": "aims",
        }
    }

    def test_correct_config_toml_only(self, tmp_path: Path, expected_config):
        config_file_path = self._make_toml_file(expected_config, tmp_path)
        actual_config = ScanConfigBuilder.build([str(config_file_path)])
        expected_config = self._make_full_expected_config(expected_config)
        self._assert_correct_config(actual_config, expected_config)

    @staticmethod
    def _assert_correct_config(actual_config, expected_config):
        assert set(actual_config.keys()) == set(expected_config.keys())
        for g in actual_config.keys():
            assert set(actual_config[g].keys()) == set(expected_config[g].keys())
            for n in actual_config[g].keys():
                if isinstance(actual_config[g][n], np.ndarray):
                    np.testing.assert_equal(actual_config[g][n], expected_config[g][n],
                                            err_msg=f'actual_config["{g}"]["{n}"] != expected_config["{g}"]["{n}"]')
                else:
                    assert actual_config[g][n] == expected_config[g][n], \
                        f'actual_config["{g}"]["{n}"] = {actual_config[g][n]}' \
                        f' != {expected_config[g][n]} = expected_config["{g}"]["{n}"]'

    def pytest_generate_tests(self, metafunc):
        if "expected_config" in metafunc.fixturenames:
            metafunc.parametrize(
                "expected_config",
                [
                    *itertools.chain(*[
                        self._gen_configs(self._build_config_example(tip_type, dft_code))
                        for tip_type, dft_code in itertools.product(self._TIP_TYPE_CHOICES,
                                                                    self._DFT_CODE_CHOICES)
                    ])
                ],
                ids=[
                    f"{tip_type}-{dft_code}-{subid}"
                    for tip_type, dft_code, subid in itertools.product(self._TIP_TYPE_CHOICES,
                                                                       self._DFT_CODE_CHOICES,
                                                                       self._gen_config_ids())
                ],
            )

    @classmethod
    def _build_config_example(cls, tip_type, dft_code):
        return {
            **cls._BASE_CONFIG_EXAMPLE,
            **cls._CONFIG_EXAMPLE_BY_TIP_TYPE[tip_type],
            **cls._CONFIG_EXAMPLE_BY_DFT_CODE[dft_code],
        }

    def _gen_configs(self, base_config: Dict[str, Any]) -> List[Dict[str, Any]]:
        return [
            {**c, **to, **so, **pbc} for c, to, so, pbc in
            itertools.product(
                [base_config],
                map(lambda to: {'tip_orb': to}, self._TIP_ORB_CHOICES),
                map(lambda so: {'sample_orbs': so}, self._SAMPLE_ORBS_CHOICES),
                map(lambda pbc: {'pbc': list(pbc)}, self._PBC_CHOICES),
            )
        ]

    def _gen_config_ids(self) -> List[str]:
        return [
            f"{to}-{so}-pbc{pbc[0]}"
            for to, so, pbc in itertools.product(
                self._TIP_ORB_CHOICES,
                self._SAMPLE_ORBS_CHOICES,
                self._PBC_CHOICES,
            )
        ]

    def _make_full_expected_config(self, base_config: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
        config = defaultdict(dict)
        for k, v in base_config.items():
            if isinstance(v, np.ndarray):
                config[self._ARGS_GROUP[k]][self._TOML_ARGS[k]] = v.tolist()
            else:
                config[self._ARGS_GROUP[k]][self._TOML_ARGS[k]] = v

        config[self._ARGS_GROUP['pbc']][self._TOML_ARGS['pbc']] = \
            list(map(float, config[self._ARGS_GROUP['pbc']][self._TOML_ARGS['pbc']]))

        if config[self._ARGS_GROUP['dft_code']][self._TOML_ARGS['dft_code']].lower() != 'cp2k':
            config[self._ARGS_GROUP['cp2k_name']][self._TOML_ARGS['cp2k_name']] = self._DEFAULT_ARG_VALUES["cp2k_name"]

        for k in self._DEFAULT_ARG_VALUES.keys():
            if self._TOML_ARGS[k] not in config[self._ARGS_GROUP[k]]:
                config[self._ARGS_GROUP[k]][self._TOML_ARGS[k]] = self._DEFAULT_ARG_VALUES[k]

        return config

    def _make_toml_file(self, config, tmp_path) -> Path:
        toml_file_path = tmp_path / 'config.toml'

        toml_config = defaultdict(dict)
        for k, v in config.items():
            if isinstance(v, Path):
                toml_config[self._ARGS_GROUP[k]][self._TOML_ARGS[k]] = str(v)
            elif isinstance(v, np.ndarray):
                toml_config[self._ARGS_GROUP[k]][self._TOML_ARGS[k]] = v.tolist()
            else:
                toml_config[self._ARGS_GROUP[k]][self._TOML_ARGS[k]] = v

        with open(toml_file_path, 'w') as f:
            toml.dump(toml_config, f)

        return toml_file_path
