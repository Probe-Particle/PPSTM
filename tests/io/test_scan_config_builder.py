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
    _ARGPARSE_ARGS = {
        "scan_type":     "--scan-type",
        "tip_type":      "--tip-type",
        "scan_window":   "--scan-window",
        "scan_dim":      "--scan-dim",
        "v":             "-v",
        "v_max":         "--v-max",
        "dv":            "-dV",
        "eta":           "--eta",
        "tip_orb":       "--tip-orb",
        "sample_orbs":   "--sample-orbs",
        "path":          "--path",
        "dft_code":      "--dft-code",
        "geometry_file": "--geometry-file",
        "pbc":           "--pbc",
        "lvs":           "--lvs",
        "cp2k_name":     "--cp2k-name",
        "q":             "-q",
        "k":             "-k",
        "data_format":   "--data-format",
        "cut_atoms":     "--cut-atoms",
        "lower_atoms":   "--lower-atoms",
        "lower_coefs":   "--lower-coefs",
        "work_function": "-wf",
        "fermi":         "--fermi",
        "cut_min":       "-cmin",
        "cut_max":       "-cmax",
        "didv_backend":  "--didv-backend",
        'plot_atoms':    "--plot-atoms",
        "save_npy":      "--save-npy",
        "save_png":      "--save-png",
        "save_wsxm":     "--save-wsxm",
        "save_xsf":      "--save-xsf",
    }

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

    _ARG_LABELS = {
        "scan_type":     "scan type",
        "tip_type":      "tip type",
        "scan_window":   "scan window",
        "scan_dim":      "scan dim",
        "v":             "V",
        "v_max":         "V max",
        "dv":            "dV",
        "eta":           "eta",
        "tip_orb":       "tip orb",
        "sample_orbs":   "sample orbs",
        "path":          "path",
        "dft_code":      "dft code",
        "geometry_file": "geometry file",
        "pbc":           "pbc",
        "lvs":           "lvs",
        "cp2k_name":     "cp2k name",
        "q":             "Q",
        "k":             "K",
        "data_format":   "data format",
        "work_function": "work function",
        "fermi":         "fermi",
        "cut_min":       "cut min",
        "cut_max":       "cut max",
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

    _MISSING_ARG_COND = defaultdict(str, {
        "scan_window": "tip_type",
        "scan_dim":    "tip_type",
        "q":           "tip_type",
        "k":           "tip_type",
        "data_format": "tip_type",
        "lvs":         "dft_code",
        "cp2k_name":   "dft_code",
    })

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

    _PBC_CHOICES = [
        [0, 0],
        [0.5, 0.5],
        [1, 1],
        [2, 2],
        [3, 3],
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

    _EXAMPLE_CONFIG_FIREBALL_DFT_FIXED_TIP = {
        "scan_type": "states",
        "tip_type": "fixed",
        "scan_window": np.asarray([[-5, -5, -5], [5, 5, 5]]),
        "scan_dim": [50, 50, 50],
        "v": -0.05,
        "v_max": 0.05,
        "dv": 0.1,
        "eta": 1e-7,
        "path": Path("./"),
        "dft_code": "fireball",
        "geometry_file": "input_plot.xyz",
        "lvs": "false",
        "work_function": 5,
        "fermi": 0,
        "cut_min": -2.5,
        "cut_max": 2.5,
    }

    _EXAMPLE_CONFIG_CP2K_DFT_FIXED_TIP = {
        "scan_type": "states",
        "tip_type": "fixed",
        "scan_window": np.asarray([[-5, -5, -5], [5, 5, 5]]),
        "scan_dim": [50, 50, 50],
        "v": -0.05,
        "v_max": 0.05,
        "dv": 0.1,
        "eta": 1e-7,
        "path": Path("./"),
        "dft_code": "cp2k",
        "geometry_file": "input_plot.xyz",
        "lvs": "false",
        "cp2k_name": "crazy_mol",
        "work_function": 5,
        "fermi": 0,
        "cut_min": -2.5,
        "cut_max": 2.5,
    }

    _EXAMPLE_CONFIG_FIREBALL_DFT_RELAXED_TIP = {
        "scan_type": "states",
        "tip_type": "relaxed",
        "v": -0.05,
        "v_max": 0.05,
        "dv": 0.1,
        "eta": 1e-7,
        "path": Path("./"),
        "dft_code": "fireball",
        "geometry_file": "input_plot.xyz",
        "lvs": "false",
        "q": 0,
        "k": 0.24,
        "data_format": "npy",
        "work_function": 5,
        "fermi": 0,
        "cut_min": -2.5,
        "cut_max": 2.5,
    }

    def test_leave_one_var_out_raise_exception(self, expected_config: Dict[str, str|List[str]], arg_name: str):
        args = self._build_args(expected_config, omitted_args=[arg_name])

        error_cond_arg = self._MISSING_ARG_COND[arg_name]

        with pytest.raises(ValueError,
                           match=rf"^The {self._ARG_LABELS[arg_name]} parameter is required"
                                 rf"{f' when {self._ARG_LABELS[error_cond_arg]} is \'{expected_config[error_cond_arg]}\'' if error_cond_arg else ''}.$"):
            ScanConfigBuilder.build(args)

    def test_correct_config_toml_only(self, tmp_path: Path, expected_config):
        config_file_path = self._make_toml_file(expected_config, tmp_path)
        actual_config = ScanConfigBuilder.build([str(config_file_path)])
        expected_config = self._make_full_expected_config(expected_config)
        self._assert_correct_config(actual_config, expected_config)

    def test_correct_config_params_only(self, expected_config: Dict[str, str | List[str]]):
        actual_config = ScanConfigBuilder.build(self._build_args(expected_config))
        expected_config = self._make_full_expected_config(expected_config)
        self._assert_correct_config(actual_config, expected_config)

    @staticmethod
    def _assert_correct_config(actual_config, expected_config):
        assert set(actual_config.keys()) == set(expected_config.keys())
        for g in actual_config.keys():
            assert set(actual_config[g].keys()) == set(expected_config[g].keys())
            for n in actual_config[g].keys():
                if isinstance(actual_config[g][n], np.ndarray):
                    np.testing.assert_equal(actual_config[g][n], expected_config[g][n])
                else:
                    assert actual_config[g][n] == expected_config[g][n]

    def pytest_generate_tests(self, metafunc):
        if "expected_config" in metafunc.fixturenames:
            if "arg_name" in metafunc.fixturenames:
                metafunc.parametrize(("expected_config", "arg_name"),
                                     [
                                         *itertools.product(
                                             self._gen_configs(self._EXAMPLE_CONFIG_FIREBALL_DFT_FIXED_TIP),
                                             [
                                                 *self._EXAMPLE_CONFIG_FIREBALL_DFT_FIXED_TIP.keys(),
                                                 'tip_orb',
                                                 'sample_orbs',
                                                 'pbc',
                                             ]),
                                         *itertools.product(
                                             self._gen_configs(self._EXAMPLE_CONFIG_FIREBALL_DFT_RELAXED_TIP),
                                             [
                                                 *self._EXAMPLE_CONFIG_FIREBALL_DFT_RELAXED_TIP.keys(),
                                                 'tip_orb',
                                                 'sample_orbs',
                                                 'pbc',
                                             ]),
                                         *itertools.product(
                                             self._gen_configs(self._EXAMPLE_CONFIG_CP2K_DFT_FIXED_TIP),
                                             [
                                                 *self._EXAMPLE_CONFIG_CP2K_DFT_FIXED_TIP.keys(),
                                                 'tip_orb',
                                                 'sample_orbs',
                                                 'pbc',
                                             ]),
                                     ])
            else:
                metafunc.parametrize("expected_config",
                                     [
                                         *self._gen_configs(self._EXAMPLE_CONFIG_FIREBALL_DFT_FIXED_TIP),
                                         *self._gen_configs(self._EXAMPLE_CONFIG_FIREBALL_DFT_RELAXED_TIP),
                                         *self._gen_configs(self._EXAMPLE_CONFIG_CP2K_DFT_FIXED_TIP),
                                     ])

    def _build_args(self, config: Dict[str, str|List[str]], omitted_args: Iterable[str]|None = None) -> List[str]:
        args = []
        if omitted_args is None:
            omitted_args = []
        for k, v in config.items():
            if k not in omitted_args:
                if k == "pbc":
                    args.extend([self._ARGPARSE_ARGS[k], str(v[0])])
                elif isinstance(v, str):
                    args.extend([self._ARGPARSE_ARGS[k], v])
                elif isinstance(v, Path):
                    args.extend([self._ARGPARSE_ARGS[k], str(v)])
                elif isinstance(v, np.ndarray):
                    args.extend([self._ARGPARSE_ARGS[k], *map(str, v.ravel().tolist())])
                elif isinstance(v, Iterable):
                    args.extend([self._ARGPARSE_ARGS[k], *map(str, v)])
                else:
                    args.extend([self._ARGPARSE_ARGS[k], str(v)])
        return args

    def _gen_configs(self, base_config: Dict[str, Any]) -> List[Dict[str, Any]]:
        return [
            {**c, **to, **so, **pbc} for c, to, so, pbc in
            itertools.product(
                [base_config],
                map(lambda to: {'tip_orb': to}, self._TIP_ORB_CHOICES),
                map(lambda so: {'sample_orbs': so}, self._SAMPLE_ORBS_CHOICES),
                map(lambda pbc: {'pbc': pbc}, self._PBC_CHOICES),
            )
        ]

    def _make_full_expected_config(self, base_config: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
        config = defaultdict(dict)
        for k, v in base_config.items():
            if isinstance(v, str) or isinstance(v, Path):
                config[self._ARGS_GROUP[k]][self._TOML_ARGS[k]] = v
            elif isinstance(v, np.ndarray):
                config[self._ARGS_GROUP[k]][self._TOML_ARGS[k]] = v.tolist()
            elif isinstance(v, Iterable):
                config[self._ARGS_GROUP[k]][self._TOML_ARGS[k]] = v
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
            if isinstance(v, str):
                toml_config[self._ARGS_GROUP[k]][self._TOML_ARGS[k]] = v
            elif isinstance(v, Path):
                toml_config[self._ARGS_GROUP[k]][self._TOML_ARGS[k]] = str(v)
            elif isinstance(v, np.ndarray):
                toml_config[self._ARGS_GROUP[k]][self._TOML_ARGS[k]] = v.tolist()
            elif isinstance(v, Iterable):
                toml_config[self._ARGS_GROUP[k]][self._TOML_ARGS[k]] = v
            else:
                toml_config[self._ARGS_GROUP[k]][self._TOML_ARGS[k]] = v

        with open(toml_file_path, 'w') as f:
            toml.dump(toml_config, f)

        return toml_file_path
