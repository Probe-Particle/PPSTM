import logging
logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)

import argparse
import os
from pathlib import Path

import numpy as np
import tomli

from ppstm import basUtils as bU
from ppstm import ReadSTM
from ppstm import STMutils
from ppstm import visualization

logger = logging.getLogger(__name__)

_MIN_TIP_SAMPLE_HEIGHT_GAP = 2.0  # Angstroms

def run_simulation(config: dict):
    # ppafm needed for relaxed tip scans and npy+xsf output
    tip_type = config['scan']['tip_type']
    npy_or_xsf_output = config['output']['NPY'] or config['output']['XSF']
    if tip_type in ['relaxed', 'r'] or npy_or_xsf_output:
        try:
            import ppafm.io as io
        except ImportError:
            raise ImportError("ppafm is required for relaxed tip scans and NPY/XSF output.")
        
    # Set tip coefficients
    tip_orb = config['scan']['tip_orb']
    tip_coeffs = STMutils.get_tip_coefficients(tip_orb)
    logger.debug(f"Set tip coefficients for {tip_orb}: {tip_coeffs}")

    # Read eigenenergies and coefficients
    eigenenergies, coefs, atoms = ReadSTM.read_dft(config)
    logger.debug(f"Read eigenenergies and coefficients.")

    # Set up atom plotting if enabled
    geom_plot = None
    if config['output']['plot_atoms']:
        try:
            geom_plot, _, _ = bU.loadAtoms(
                os.path.join(
                    config['input']['path'],
                    'input_plot.xyz'
                )
            ) 
        except FileNotFoundError:
            logger.warning("Atom plotting disabled due to missing `input_plot.xyz` file.")

    # Get tip positions
    (
        tip_r,
        tip_r0,
        lvec,
        extent,
        atomic_head_or_info
    ) = STMutils.get_tip_positions(config)
    logger.debug(f"Tip positions read for a {config['scan']['tip_type']} scan.")

    if not _is_tip_above_sample(tip_positions=tip_r,
                                sample_atom_positions=atoms,
                                min_gap=_MIN_TIP_SAMPLE_HEIGHT_GAP):
        logger.warning(f"Detected minimum tip height below the maximum sample height "
                       f"within {_MIN_TIP_SAMPLE_HEIGHT_GAP} Å")

    # Run STM scan
    current, didv = STMutils.run_stm_scan(
        config,
        eigenenergies,
        coefs,
        tip_r,
        atoms
    )
    logger.debug(f"STM scan complete for scan type {config['scan']['scan_type']}.")

    # Get voltages and names
    voltages, names = visualization.get_voltages_and_names(config, eigenenergies)

    # Output results
    if config['output']['PNG']:
        atom_size = config['output'].get('atom_size', 0.15)
        visualization.plot_png(config, current, didv, voltages, names, lvec, extent, geom_plot, atom_size)
        logger.debug(f"PNG output complete.")
    if config['output']['WSxM']:
        visualization.plot_wsxm(config, current, didv, voltages, names, tip_r0)
        logger.debug(f"WSXM output complete.")
    if config['output']['XSF']:
        visualization.save_xsf(config, current, didv, voltages, names, geom_plot, lvec)
        logger.debug(f"XSF output complete.")
    if config['output']['NPY']:
        visualization.save_npz(config, current, didv, voltages, names, lvec, atomic_head_or_info)
        logger.debug(f"NPY output complete.")

    logger.debug(f"Output finished, exiting.")

def _is_tip_above_sample(tip_positions: np.ndarray, sample_atom_positions: np.ndarray, min_gap: float) -> bool:
    """Checks if the lowest tip point is safely above the highest sample point.

    Precondition: Both input arrays must be non-empty.
    """
    lowest_tip_height = np.min(tip_positions[...,2])
    highest_sample_height = np.max(sample_atom_positions[:,2])
    return bool(lowest_tip_height >= highest_sample_height + min_gap)

def _existing_toml_file(value: str) -> Path:
    path = Path(value)

    if not path.is_file():
        raise argparse.ArgumentTypeError(f"file does not exist: {path}")

    if path.suffix.lower() != ".toml":
        raise argparse.ArgumentTypeError(f"expected a .toml file: {path}")

    return path

def main():
    parser = argparse.ArgumentParser(
        description="Execute PP-STM simulation scan",
    )
    parser.add_argument(
        'config_file',
        type=_existing_toml_file,
        help="TOML configuration file with PP-STM simulation parameters",
    )
    args = parser.parse_args()

    # Get config file from command line
    config_file = args.config_file

    # Load config file
    with open(config_file, 'rb') as f:
        config = tomli.load(f)

    logger.debug(f"Loaded config from {config_file}")
    logger.debug(f"Config: {config}")

    run_simulation(config)

if __name__=='__main__':
    main()
