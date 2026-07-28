import logging
import os
from typing import Iterable

import numpy as np

from pyPPSTM import ReadSTM
from pyPPSTM import STMutils
from pyPPSTM import basUtils as bU
from pyPPSTM import visualization
from pyPPSTM.scan_config_builder import ScanConfigBuilder

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
    print(f"Set tip coefficients for {tip_orb}: {tip_coeffs}")

    # Read eigenenergies and coefficients
    eigenenergies, coefs, atoms = ReadSTM.read_dft(config)
    print(f"Read eigenenergies and coefficients.")

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
            print("WARNING: Atom plotting disabled due to missing input_plot.xyz file.")

    # Get tip positions
    (
        tip_r,
        tip_r0,
        lvec,
        extent,
        atomic_head_or_info
    ) = STMutils.get_tip_positions(config)
    print(f"Tip positions read for a {config['scan']['tip_type']} scan.")

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
    print(f"STM scan complete for scan type {config['scan']['scan_type']}.")

    # Get voltages and names
    voltages, names = visualization.get_voltages_and_names(config, eigenenergies)

    # Output results
    if config['output']['PNG']:
        atom_size = config['output'].get('atom_size', 0.15)
        visualization.plot_png(config, current, didv, voltages, names, lvec, extent, geom_plot, atom_size)
        print(f"PNG output complete.")
    if config['output']['WSxM']:
        visualization.plot_wsxm(config, current, didv, voltages, names, tip_r0)
        print(f"WSXM output complete.")
    if config['output']['XSF']:
        visualization.save_xsf(config, current, didv, voltages, names, geom_plot, lvec)
        print(f"XSF output complete.")
    if config['output']['NPY']:
        visualization.save_npz(config, current, didv, voltages, names, lvec, atomic_head_or_info)
        print(f"NPY output complete.")

    print(f"Output finished, exiting.")

def _is_tip_above_sample(tip_positions: np.ndarray, sample_atom_positions: np.ndarray, min_gap: float) -> bool:
    """Checks if the lowest tip point is safely above the highest sample point.

    Precondition: Both input arrays must be non-empty.
    """
    lowest_tip_height = np.min(tip_positions[...,2])
    highest_sample_height = np.max(sample_atom_positions[:,2])
    return bool(lowest_tip_height >= highest_sample_height + min_gap)


def main(args: Iterable[str]|None = None):
    config = ScanConfigBuilder.build(args)

    run_simulation(config)


if __name__=='__main__':
    logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)
    main()
