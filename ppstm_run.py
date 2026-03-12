import os
import sys
import tomli

from pyPPSTM import basUtils as bU
from pyPPSTM import ReadSTM
from pyPPSTM import STMutils
from pyPPSTM import visualization

def main(config: dict):
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
    if config['output']['plot_atoms']:
        try:
            geom_plot, _, _ = bU.loadAtoms(
                os.path.join(
                    config['input']['path'],
                    'input_plot.xyz'
                )
            ) 
        except FileNotFoundError:
            geom_plot = None
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
        visualization.plot_png(config, current, didv, voltages, names, lvec, extent, geom_plot)
        print(f"PNG output complete.")
    if config['output']['WSxM']:
        visualization.plot_wsxm(config, current, didv, voltages, names, tip_r0)
        print(f"WSXM output complete.")
    if config['output']['XSF']:
        visualization.save_xsf(config, current, didv, voltages, names, geom_plot, lvec)
        print(f"XSF output complete.")
    if config['output']['NPY']:
        visualization.save_npy(config, current, didv, voltages, names, lvec, atomic_head_or_info)
        print(f"NPY output complete.")

    print(f"Output finished, exiting.")

if __name__=='__main__':
    # Get config file from command line
    config_file = sys.argv[1]

    # Load config file
    with open(config_file, 'rb') as f:
        config = tomli.load(f)
    
    print(f"Loaded config from {config_file}")
    print(f"Config: {config}")

    main(config)
