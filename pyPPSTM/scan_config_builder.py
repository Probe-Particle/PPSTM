import argparse
import logging
from pathlib import Path
from typing import Iterable

import numpy as np
import tomli

from pyPPSTM import DidvBackend


class ScanConfigBuilder:
    _DEFAULT_CP2K_NAME = ''
    _DEFAULT_CUT_ATOMS = -1
    _DEFAULT_LOWER_COEFS = []
    _DEFAULT_LOWER_ATOMS = []
    _DEFAULT_DIDV_BACKEND = DidvBackend.CPP

    @classmethod
    def build(cls, args: Iterable[str] | None = None):
        args = cls._parse_args(args)

        if args.config_file:
            # Load config file
            with open(args.config_file, 'rb') as f:
                config = tomli.load(f)
        else:
            config = dict()

        config = cls._check_configuration(config, args)

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
        parser.add_argument(
            '-s', '--scan-type',
            type=str,
            choices=[
                'didv', 'dIdV', 'didv-single',
                'states', 'STATES',
                'STM', 'STM-single',
                'v-scan', 'V-scan', 'Voltage-scan',
            ],
        )
        parser.add_argument(
            '-t', '--tip-type',
            type=str,
            choices=[
                'relaxed', 'r',
                'fixed', 'f',
            ],
        )
        parser.add_argument(
            '-w', '--scan-window',
            action="extend",
            nargs="+",
            type=float,
        )
        parser.add_argument(
            '-d', '--scan-dim',
            action="extend",
            nargs="+",
            type=float,
        )
        parser.add_argument(
            '-v', '-V',
            type=float,
        )
        parser.add_argument(
            '-vM', '-VM', '--v-max',
            type=float,
        )
        parser.add_argument(
            '-dv', '-dV',
            type=float,
        )
        parser.add_argument(
            '--eta',
            type=float,
        )
        parser.add_argument(
            '-to', '--tip-orb',
            type=str,
            choices=['s', 'pxy', 'spxy', '5spxy', '10spxy', 'CO', 'pz', 'dz2', 'dxzyz'],
        )
        parser.add_argument(
            '-so', '--sample-orbs',
            type=str,
            choices=['sp', 'spd'],
        )
        parser.add_argument(
            '-p', '--path',
            type=Path,
        )
        parser.add_argument(
            '--dft-code',
            type=str,
            choices=[
                'fireball', 'Fireball', 'FIREBALL',
                'aims', 'AIMS', 'FHI-AIMS',
                'cp2k', 'CP2K',
                'gpaw', 'GPAW'
            ],
        )
        parser.add_argument(
            '-gf', '--geometry-file',
            type=Path,
        )
        parser.add_argument(
            '-pbc', '--pbc',
            type=float,
            choices=[0, 0.5, 1, 2, 3],
        )
        parser.add_argument(
            '-lvs', '--lvs',
        )
        parser.add_argument(
            '--spin',
            type=str,
        )
        parser.add_argument(
            '--cp2k-name',
            type=str,
        )
        parser.add_argument(
            '-q', '-Q',
            type=float,
        )
        parser.add_argument(
            '-k', '-K',
            type=float,
        )
        parser.add_argument(
            '-df', '--data-format',
            type=str,
            choices=[
                'xsf', 'XSF',
                'npy', 'NPY',
            ]
        )
        parser.add_argument(
            '-ca', '--cut-atoms',
            type=cls._check_cut_atoms,
            default=cls._DEFAULT_CUT_ATOMS,
        )
        parser.add_argument(
            '-la', '--lower-atoms',
            action="extend",
            nargs="+",
            type=int,
            default=cls._DEFAULT_LOWER_ATOMS,
        )
        parser.add_argument(
            '-lc', '--lower-coefs',
            action="extend",
            nargs="+",
            type=float,
            default=cls._DEFAULT_LOWER_COEFS,
        )
        parser.add_argument(
            '-wf', '--work-function',
            type=float,
        )
        parser.add_argument(
            '-wfd', '--work-function-decay',
            type=float,
        )
        parser.add_argument(
            '-f', '--fermi',
            type=float,
        )
        parser.add_argument(
            '-cm', '-cmin', '--cut-min',
            type=float,
        )
        parser.add_argument(
            '-cM', '-cmax', '--cut-max',
            type=float,
        )
        parser.add_argument(
            '-be', '--didv-backend',
            type=DidvBackend,
            choices=list(DidvBackend),
            default=cls._DEFAULT_DIDV_BACKEND,
        )
        parser.add_argument(
            '--plot-atoms',
            action="store_true",
        )
        parser.add_argument(
            '-npy', '--save-npy',
            action="store_true",
        )
        parser.add_argument(
            '-xsf', '--save-xsf',
            action="store_true",
        )
        parser.add_argument(
            '-png', '--save-png',
            action="store_true",
        )
        parser.add_argument(
            '-wsxm', '--save-wsxm',
            action="store_true",
        )
        args = parser.parse_args(args)

        if args.config_file:
            cls._existing_toml_file(args.config_file)

        return args

    @staticmethod
    def _check_cut_atoms(value: str) -> int:
        value = int(value)

        if value == 0 or value < -1:
            raise argparse.ArgumentTypeError(f"cut atoms should be greater than 0 or -1")

        return value

    @staticmethod
    def _existing_toml_file(value: str) -> Path:
        path = Path(value)

        if not path.is_file():
            raise argparse.ArgumentTypeError(f"file does not exist: {path}")

        if path.suffix.lower() != ".toml":
            raise argparse.ArgumentTypeError(f"expected a .toml file: {path}")

        return path

    @classmethod
    def _check_configuration(cls, config, args):
        for g in ['scan', 'input', 'advanced', 'output']:
            if g not in config:
                config[g] = {}

        if args.scan_type is not None:
            config['scan']['scan_type'] = args.scan_type
        elif 'scan_type' not in config['scan']:
            raise ValueError("The scan type parameter is required.")

        if args.tip_type is not None:
            config['scan']['tip_type'] = args.tip_type
        elif 'tip_type' not in config['scan']:
            raise ValueError("The tip type parameter is required.")

        if args.scan_window is not None and len(args.scan_window) > 0:
            if len(args.scan_window) < 6:
                raise argparse.ArgumentTypeError(f"expected six numbers to define the scan window, "
                                                 f"got only {len(args.scan_window)} values.")
            config['scan']['scan_window'] = args.scan_window
        if 'scan_window' in config['scan']:
            config['scan']['scan_window'] = np.asarray(config['scan']['scan_window']).reshape(2, 3)
        elif config['scan']['tip_type'] in ['fixed', 'f']:
            raise ValueError(f"The scan window parameter is required "
                             f"when tip type is '{config['scan']['tip_type']}'.")

        if args.scan_dim is not None and len(args.scan_dim) > 0:
            if len(args.scan_dim) < 3:
                raise argparse.ArgumentTypeError(f"expected three numbers to define the scan dim, "
                                                 f"got only {len(args.scan_dim)} values.")
            config['scan']['scan_dim'] = args.scan_dim
        elif 'scan_dim' not in config['scan'] and config['scan']['tip_type'] in ['fixed', 'f']:
            raise ValueError(f"The scan dim parameter is required "
                             f"when tip type is '{config['scan']['tip_type']}'.")

        if args.v is not None:
            config['scan']['V'] = args.v
        elif 'V' not in config['scan']:
            raise ValueError("The V parameter is required.")

        if args.v_max is not None:
            config['scan']['V_max'] = args.v_max
        elif 'V_max' not in config['scan']:
            raise ValueError("The V max parameter is required.")

        if args.dv is not None:
            config['scan']['dV'] = args.dv
        elif 'dV' not in config['scan']:
            raise ValueError("The dV parameter is required.")

        if args.eta is not None:
            config['scan']['eta'] = args.eta
        elif 'eta' not in config['scan']:
            raise ValueError("The eta parameter is required.")

        if args.tip_orb is not None:
            config['scan']['tip_orb'] = args.tip_orb
        elif 'tip_orb' not in config['scan']:
            raise ValueError("The tip orb parameter is required.")

        if args.sample_orbs is not None:
            config['scan']['sample_orbs'] = args.sample_orbs
        elif 'sample_orbs' not in config['scan']:
            raise ValueError("The sample orbs parameter is required.")

        if args.path is not None:
            config['input']['path'] = args.path
        elif 'path' in config['input']:
            config['input']['path'] = Path(config['input']['path'])
        else:
            raise ValueError("The path parameter is required.")

        if args.dft_code is not None:
            config['input']['dft_code'] = args.dft_code
        elif 'dft_code' not in config['input']:
            raise ValueError("The dft code parameter is required.")

        if args.geometry_file is not None:
            config['input']['geometry_file'] = str(args.geometry_file)
        elif 'geometry_file' not in config['input']:
            raise ValueError("The geometry file parameter is required.")

        if args.pbc is not None:
            config['input']['pbc'] = [args.pbc, ] * 2
        elif 'pbc' not in config['input']:
            raise ValueError("The pbc parameter is required.")

        if args.lvs is not None:
            config['input']['lvs'] = args.lvs
        elif 'lvs' not in config['input'] and config['input']['dft_code'].lower() in ('fireball', 'cp2k'):
            raise ValueError("The lvs parameter is required "
                             f"when dft code is '{config['input']['dft_code']}'.")

        if args.spin is not None or 'spin' not in config['input']:
            config['input']['spin'] = args.spin

        if config['input']['dft_code'].lower() == 'cp2k':
            if args.cp2k_name is not None:
                config['input']['cp2k_name'] = args.cp2k_name
            elif 'cp2k_name' not in config['input']:
                raise ValueError("The cp2k name parameter is required "
                                 f"when dft code is '{config['input']['dft_code']}'.")
        else:
            config['input']['cp2k_name'] = cls._DEFAULT_CP2K_NAME

        if args.q is not None:
            config['input']['Q'] = args.q
        elif 'Q' not in config['input'] and config['scan']['tip_type'] in ['relaxed', 'r']:
            raise ValueError(f"The Q parameter is required "
                             f"when tip type is '{config['scan']['tip_type']}'.")

        if args.k is not None:
            config['input']['K'] = args.k
        elif 'K' not in config['input'] and config['scan']['tip_type'] in ['relaxed', 'r']:
            raise ValueError(f"The K parameter is required "
                             f"when tip type is '{config['scan']['tip_type']}'.")

        if args.data_format is not None:
            config['input']['data_format'] = args.data_format
        elif 'data_format' not in config['input'] and config['scan']['tip_type'] in ['relaxed', 'r']:
            raise ValueError(f"The data format parameter is required "
                             f"when tip type is '{config['scan']['tip_type']}'.")

        if 'cut_atoms' not in config['advanced']:
            config['advanced']['cut_atoms'] = args.cut_atoms

        if args.lower_atoms is not None:
            config['advanced']['lower_atoms'] = args.lower_atoms

        if args.lower_coefs is not None:
            config['advanced']['lower_coefs'] = args.lower_coefs

        if len(config['advanced']['lower_atoms']) != len(config['advanced']['lower_coefs']):
            raise ValueError(
                "lower atoms and lower coefs must have the same length: "
                f"got {len(config['advanced']['lower_atoms'])} lower atoms value(s) "
                f"and {len(config['advanced']['lower_coefs'])} lower coefs value(s). "
                "Each lower atoms entry must have a corresponding lower coefs entry."
            )

        if args.work_function is not None:
            config['advanced']['work_function'] = args.work_function
        elif 'work_function' not in config['advanced']:
            raise ValueError("The work function parameter is required.")

        if args.work_function_decay is not None:
            config['advanced']['work_function_decay'] = args.work_function_decay
        elif 'work_function_decay' not in config['advanced'] and \
                config['scan']['scan_type'] in ['STM', 'STM-single', 'v-scan', 'V-scan', 'Voltage-scan']:
            raise ValueError(f"The work function decay parameter is required "
                             f"when scan type is '{config['scan']['scan_type']}'.")

        if args.fermi is not None:
            config['advanced']['fermi'] = args.fermi
        elif 'fermi' not in config['advanced']:
            raise ValueError("The fermi parameter is required.")

        if args.cut_min is not None:
            config['advanced']['cut_min'] = args.cut_min
        elif 'cut_min' not in config['advanced']:
            raise ValueError("The cut min parameter is required.")

        if args.cut_max is not None:
            config['advanced']['cut_max'] = args.cut_max
        elif 'cut_max' not in config['advanced']:
            raise ValueError("The cut max parameter is required.")

        if args.didv_backend is not None:
            config['advanced']['didv_backend'] = args.didv_backend

        if args.plot_atoms is not None or 'plot_atoms' not in config['output']:
            config['output']['plot_atoms'] = args.plot_atoms

        if args.save_npy is not None or 'NPY' not in config['output']:
            config['output']['NPY'] = args.save_npy

        if args.save_xsf is not None or 'XSF' not in config['output']:
            config['output']['XSF'] = args.save_xsf

        if args.save_png is not None or 'PNG' not in config['output']:
            config['output']['PNG'] = args.save_png

        if args.save_wsxm is not None or 'WSxM' not in config['output']:
            config['output']['WSxM'] = args.save_wsxm

        return config
