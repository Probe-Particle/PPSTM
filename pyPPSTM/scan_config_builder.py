import argparse
import logging
from pathlib import Path
from typing import Iterable, Dict, Any

import numpy as np
import tomli

from pyPPSTM import DidvBackend


class ScanConfigBuilder:
    _DEFAULT_PLOT_ATOMS = _DEFAULT_NPY = _DEFAULT_XSF = _DEFAULT_PNG = _DEFAULT_WSXM = False
    _DEFAULT_SPIN = None
    _DEFAULT_CP2K_NAME = ''
    _DEFAULT_CUT_ATOMS = -1
    _DEFAULT_LOWER_COEFS = []
    _DEFAULT_LOWER_ATOMS = []
    _DEFAULT_DIDV_BACKEND = DidvBackend.CPP

    @classmethod
    def build(cls, args: Iterable[str] | None = None) -> Dict[str, Dict[str, Any]]:
        args = cls._parse_args(args)

        if args.config_file:
            # Load config file
            with open(args.config_file, 'rb') as f:
                config = tomli.load(f)
        else:
            config = dict()

        config = cls._check_configuration(config)

        logging.info(f"Config: {config}")

        return config

    @classmethod
    def _parse_args(cls, args: Iterable[str] | None = None):
        parser = argparse.ArgumentParser(
            description="Execute PP-STM simulation scan",
        )
        parser.add_argument(
            'config_file',
            type=Path,
            nargs='?',
            help="TOML configuration file with PP-STM simulation parameters",
        )
        args = parser.parse_args(args)

        cls._existing_toml_file(args.config_file)

        return args

    @staticmethod
    def _existing_toml_file(value: str) -> Path:
        path = Path(value)

        if not path.is_file():
            raise argparse.ArgumentTypeError(f"file does not exist: {path}")

        if path.suffix.lower() != ".toml":
            raise argparse.ArgumentTypeError(f"expected a .toml file: {path}")

        return path

    @classmethod
    def _check_configuration(cls, config) -> Dict[str, Dict[str, Any]]:
        for g in ['scan', 'input', 'advanced', 'output']:
            if g not in config:
                config[g] = {}

        if 'scan_type' not in config['scan']:
            raise ValueError("The scan type parameter is required.")

        if 'tip_type' not in config['scan']:
            raise ValueError("The tip type parameter is required.")

        if 'scan_window' in config['scan']:
            config['scan']['scan_window'] = np.asarray(config['scan']['scan_window']).reshape(2, 3)
        elif config['scan']['tip_type'] in ['fixed', 'f']:
            raise ValueError(f"The scan window parameter is required "
                             f"when tip type is '{config['scan']['tip_type']}'.")

        if 'scan_dim' not in config['scan'] and config['scan']['tip_type'] in ['fixed', 'f']:
            raise ValueError(f"The scan dim parameter is required "
                             f"when tip type is '{config['scan']['tip_type']}'.")

        if 'V' not in config['scan']:
            raise ValueError("The V parameter is required.")

        if 'V_max' not in config['scan']:
            raise ValueError("The V max parameter is required.")

        if 'dV' not in config['scan']:
            raise ValueError("The dV parameter is required.")

        if 'eta' not in config['scan']:
            raise ValueError("The eta parameter is required.")

        if 'tip_orb' not in config['scan']:
            raise ValueError("The tip orb parameter is required.")

        if 'sample_orbs' not in config['scan']:
            raise ValueError("The sample orbs parameter is required.")

        if 'path' in config['input']:
            config['input']['path'] = Path(config['input']['path'])
        else:
            raise ValueError("The path parameter is required.")

        if 'dft_code' not in config['input']:
            raise ValueError("The dft code parameter is required.")

        if 'geometry_file' not in config['input']:
            raise ValueError("The geometry file parameter is required.")

        if 'pbc' not in config['input']:
            raise ValueError("The pbc parameter is required.")

        if 'lvs' not in config['input'] and config['input']['dft_code'].lower() in ('fireball', 'cp2k'):
            raise ValueError("The lvs parameter is required "
                             f"when dft code is '{config['input']['dft_code']}'.")

        if 'spin' not in config['input']:
            config['input']['spin'] = cls._DEFAULT_SPIN

        if config['input']['dft_code'].lower() == 'cp2k':
            if 'cp2k_name' not in config['input']:
                raise ValueError("The cp2k name parameter is required "
                                 f"when dft code is '{config['input']['dft_code']}'.")
        else:
            config['input']['cp2k_name'] = cls._DEFAULT_CP2K_NAME

        if 'Q' not in config['input'] and config['scan']['tip_type'] in ['relaxed', 'r']:
            raise ValueError(f"The Q parameter is required "
                             f"when tip type is '{config['scan']['tip_type']}'.")

        if 'K' not in config['input'] and config['scan']['tip_type'] in ['relaxed', 'r']:
            raise ValueError(f"The K parameter is required "
                             f"when tip type is '{config['scan']['tip_type']}'.")

        if 'data_format' not in config['input'] and config['scan']['tip_type'] in ['relaxed', 'r']:
            raise ValueError(f"The data format parameter is required "
                             f"when tip type is '{config['scan']['tip_type']}'.")

        if 'cut_atoms' not in config['advanced']:
            config['advanced']['cut_atoms'] = cls._DEFAULT_CUT_ATOMS

        if 'lower_atoms' not in config['advanced']:
            config['advanced']['lower_atoms'] = cls._DEFAULT_LOWER_ATOMS

        if 'lower_coefs' not in config['advanced']:
            config['advanced']['lower_coefs'] = cls._DEFAULT_LOWER_COEFS

        if len(config['advanced']['lower_atoms']) != len(config['advanced']['lower_coefs']):
            raise ValueError(
                "lower atoms and lower coefs must have the same length: "
                f"got {len(config['advanced']['lower_atoms'])} lower atoms value(s) "
                f"and {len(config['advanced']['lower_coefs'])} lower coefs value(s). "
                "Each lower atoms entry must have a corresponding lower coefs entry."
            )

        if 'work_function' not in config['advanced']:
            raise ValueError("The work function parameter is required.")

        if 'work_function_decay' not in config['advanced'] and \
                config['scan']['scan_type'] in ['STM', 'STM-single', 'v-scan', 'V-scan', 'Voltage-scan']:
            raise ValueError(f"The work function decay parameter is required "
                             f"when scan type is '{config['scan']['scan_type']}'.")

        if 'fermi' not in config['advanced']:
            raise ValueError("The fermi parameter is required.")

        if 'cut_min' not in config['advanced']:
            raise ValueError("The cut min parameter is required.")

        if 'cut_max' not in config['advanced']:
            raise ValueError("The cut max parameter is required.")

        if 'didv_backend' not in config['advanced']:
            config['advanced']['didv_backend'] = cls._DEFAULT_DIDV_BACKEND

        if 'plot_atoms' not in config['output']:
            config['output']['plot_atoms'] = cls._DEFAULT_PLOT_ATOMS

        if 'NPY' not in config['output']:
            config['output']['NPY'] = cls._DEFAULT_NPY

        if 'XSF' not in config['output']:
            config['output']['XSF'] = cls._DEFAULT_XSF

        if 'PNG' not in config['output']:
            config['output']['PNG'] = cls._DEFAULT_PNG

        if 'WSxM' not in config['output']:
            config['output']['WSxM'] = cls._DEFAULT_WSXM

        return config
