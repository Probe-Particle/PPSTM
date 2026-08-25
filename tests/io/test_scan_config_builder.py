import itertools
from collections import defaultdict
from pathlib import Path
from typing import Iterable, List, Dict, Any

import numpy as np
import pytest
import toml

from pyPPSTM import DidvBackend
from pyPPSTM.scan_config_builder import ScanConfigBuilder


class TestScanConfigBuilder:
    _CHOICES = {
        "tip_type": [
            "fixed",
            "relaxed",
        ],
        "tip_orb": [
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
        ],
        "sample_orbs": [
            "sp",
            "spd",
        ],
        "dft_code": [
            "fireball",
            "cp2k",
            "gpaw",
            "aims",
        ],
        "pbc": [
            [0, 0],
            [0.5, 0.5],
            [1, 1],
            [2, 2],
            [3, 3],
        ],
    }

    _TOML_ARG_BY_ID = {
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

    _ARG_LABEL_BY_ID = {
        "scan_type": "scan type",
        "tip_type": "tip type",
        "scan_window": "scan window",
        "scan_dim": "scan dim",
        "v": "V",
        "v_max": "V max",
        "dv": "dV",
        "eta": "eta",
        "tip_orb": "tip orb",
        "sample_orbs": "sample orbs",
        "path": "path",
        "dft_code": "dft code",
        "geometry_file": "geometry file",
        "pbc": "pbc",
        "lvs": "lvs",
        "cp2k_name": "cp2k name",
        "q": "Q",
        "k": "K",
        "data_format": "data format",
        "work_function": "work function",
        "fermi": "fermi",
        "cut_min": "cut min",
        "cut_max": "cut max",
    }

    _ARG_GROUP_BY_ARG_ID = {
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

    _REQUIRED_ARG_COND = defaultdict(str, {
        "scan_window": "tip_type",
        "scan_dim": "tip_type",
        "q": "tip_type",
        "k": "tip_type",
        "data_format": "tip_type",
        "lvs": "dft_code",
        "cp2k_name": "dft_code",
    })

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

    def test_correct_config(self, tmp_path: Path, expected_config: Dict[str, str|List[str]]):
        config_file_path = self._make_toml_file(expected_config, tmp_path)
        actual_config = ScanConfigBuilder.build([str(config_file_path)])
        expected_config = self._make_full_expected_config(expected_config)
        self._assert_correct_config(actual_config, expected_config)

    def test_leave_one_var_out_raise_exception(self,
                                               tmp_path: Path,
                                               expected_config: Dict[str, str|List[str]],
                                               arg_name: str):
        config_file_path = self._make_toml_file(expected_config, tmp_path, omitted_args=[arg_name])

        with pytest.raises(ValueError,
                           match=rf"^The {self._ARG_LABEL_BY_ID[arg_name]} parameter is required"
                                 rf"{f' when {self._ARG_LABEL_BY_ID[required_arg_cond]}'
                                     f' is \'{expected_config[required_arg_cond]}\'' 
                                 if (required_arg_cond := self._REQUIRED_ARG_COND[arg_name]) else ''}.$"
                           ):
            ScanConfigBuilder.build([str(config_file_path)])

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
            if "arg_name" in metafunc.fixturenames:
                metafunc.parametrize(
                    ("expected_config", "arg_name"),
                    [
                        *itertools.chain(*[
                            itertools.product(
                                self._gen_configs(self._build_config_example(tip_type, dft_code)),
                                self._gen_arg_names(self._build_config_example(tip_type, dft_code)),
                            )
                            for tip_type, dft_code in itertools.product(self._CHOICES["tip_type"],
                                                                        self._CHOICES["dft_code"])
                        ])
                    ],
                    ids=itertools.chain(*[
                        [
                            f"{tip_type}-{dft_code}-{subid}-without-{arg_name}"
                            for arg_name in self._gen_arg_names(self._build_config_example(tip_type, dft_code))
                        ]
                        for tip_type, dft_code, subid in itertools.product(self._CHOICES["tip_type"],
                                                                           self._CHOICES["dft_code"],
                                                                           self._gen_config_ids())
                    ])
                )
            else:
                metafunc.parametrize(
                    "expected_config",
                    [
                        *itertools.chain(*[
                            self._gen_configs(self._build_config_example(tip_type, dft_code))
                            for tip_type, dft_code in itertools.product(self._CHOICES["tip_type"],
                                                                        self._CHOICES["dft_code"])
                        ])
                    ],
                    ids=[
                        f"{tip_type}-{dft_code}-{subid}"
                        for tip_type, dft_code, subid in itertools.product(self._CHOICES["tip_type"],
                                                                           self._CHOICES["dft_code"],
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
                map(lambda to: {'tip_orb': to}, self._CHOICES["tip_orb"]),
                map(lambda so: {'sample_orbs': so}, self._CHOICES["sample_orbs"]),
                map(lambda pbc: {'pbc': pbc}, self._CHOICES["pbc"]),
            )
        ]

    def _gen_config_ids(self) -> List[str]:
        return [
            f"{to}-{so}-pbc{pbc[0]}"
            for to, so, pbc in itertools.product(
                self._CHOICES["tip_orb"],
                self._CHOICES["sample_orbs"],
                self._CHOICES["pbc"],
            )
        ]

    def _gen_arg_names(self, base_config: Dict[str, Any]) -> List[str]:
        return [
            *base_config.keys(),
            'tip_orb',
            'sample_orbs',
            'pbc',
        ]

    def _make_full_expected_config(self, base_config: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
        config = defaultdict(dict)
        for k, v in base_config.items():
            if isinstance(v, np.ndarray):
                config[self._ARG_GROUP_BY_ARG_ID[k]][self._TOML_ARG_BY_ID[k]] = v.tolist()
            else:
                config[self._ARG_GROUP_BY_ARG_ID[k]][self._TOML_ARG_BY_ID[k]] = v

        config[self._ARG_GROUP_BY_ARG_ID['pbc']][self._TOML_ARG_BY_ID['pbc']] = \
            list(map(float, config[self._ARG_GROUP_BY_ARG_ID['pbc']][self._TOML_ARG_BY_ID['pbc']]))

        if config[self._ARG_GROUP_BY_ARG_ID['dft_code']][self._TOML_ARG_BY_ID['dft_code']].lower() != 'cp2k':
            config[self._ARG_GROUP_BY_ARG_ID['cp2k_name']][self._TOML_ARG_BY_ID['cp2k_name']] = self._DEFAULT_ARG_VALUES["cp2k_name"]

        for k in self._DEFAULT_ARG_VALUES.keys():
            if self._TOML_ARG_BY_ID[k] not in config[self._ARG_GROUP_BY_ARG_ID[k]]:
                config[self._ARG_GROUP_BY_ARG_ID[k]][self._TOML_ARG_BY_ID[k]] = self._DEFAULT_ARG_VALUES[k]

        return config

    def _make_toml_file(self, config: Dict[str, Any], tmp_path: Path, omitted_args: Iterable[str]|None = None) -> Path:
        if omitted_args is None:
            omitted_args = ()
        toml_file_path = tmp_path / 'config.toml'
        toml_config = defaultdict(dict)
        for k, v in config.items():
            if k not in omitted_args:
                if isinstance(v, Path):
                    toml_config[self._ARG_GROUP_BY_ARG_ID[k]][self._TOML_ARG_BY_ID[k]] = str(v)
                elif isinstance(v, np.ndarray):
                    toml_config[self._ARG_GROUP_BY_ARG_ID[k]][self._TOML_ARG_BY_ID[k]] = v.tolist()
                else:
                    toml_config[self._ARG_GROUP_BY_ARG_ID[k]][self._TOML_ARG_BY_ID[k]] = v

        with open(toml_file_path, 'w') as f:
            toml.dump(toml_config, f)

        return toml_file_path
