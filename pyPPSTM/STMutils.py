import os
import numpy as np

from typing import Dict, Tuple

import pyPPSTM.ReadSTM as RS
import pyPPSTM.ProbeSTM as PS
import pyPPSTM.basUtils as BU

def get_tip_coefficients(tip_orb: str) -> Dict:
    """
    Returns list of coefficients for tip orbitals.
    's' -- s,
    'pxy' -- px & py,
    'spxy' -- 50% s & 50% pxy,
    '5spxy' -- 5% s & 95% pxy,
    '10spxy' -- 10% s & 90% pxy,
    'CO' -- 13% s & 87% pxy (PRL 119, 166001 (2017)),
    'pz' -- pz,
    For sample_orbs = 'sp' , possible 'dz2' and 'dxzyz' -- dxz & dyz

    Args:
        tip_orb (str): The tip orbital.

    Returns:
        coeffs (Dict): The list of coefficients for the tip orbital.    
    """
    tc = {
        "s": 0.0,
        "px": 0.0,
        "py": 0.0,
        "pz": 0.0,
        "dz2": 0.0,
        "dxz": 0.0,
        "dyz": 0.0
    }
    # [s, px, py, pz, dz2, dxz, dyz ]
    if (tip_orb == 's'):
        tc['s'] = 1.0
    elif (tip_orb == 'pxy'):
        tc['px'] = 0.5
        tc['py'] = 0.5
    elif (tip_orb == 'spxy'):
        tc['s'] = 0.5
        tc['px'] = 0.25
        tc['py'] = 0.25
    elif (tip_orb == '5spxy'):
        tc['s'] = 0.05
        tc['px'] = 0.475
        tc['py'] = 0.475
    elif (tip_orb == '10spxy'):
        tc['s'] = 0.10
        tc['px'] = 0.45
        tc['py'] = 0.45
    elif (tip_orb == 'CO'):
        tc['s'] = 0.15
        tc['px'] = 0.5
        tc['py'] = 0.5
    elif (tip_orb == 'pz'):
        tc['pz'] = 1.0
    elif (tip_orb == 'dz2'):
        tc['dz2'] = 1.0
    elif (tip_orb == 'dxzyz'):
        tc['dxz'] = 0.5
        tc['dyz'] = 0.5
    else:
        raise ValueError(f"Invalid tip orbital: {tip_orb}")
    return tc

def get_tip_positions(config: dict):
    """
    Get the tip positions for relaxed or fixed scan.

    Args:
        config (dict): The configuration dictionary.

    Returns:
        tip_r (np.ndarray): The tip positions.
        tip_r0 (np.ndarray): The initial tip positions.
        lvec (np.ndarray): The lattice vectors.
        extent (tuple): The extent of the scan.
        atomic_info_or_head (tuple): The atomic information or head.
    """
    tip_type = config['scan']['tip_type']

    if ((tip_type =='relaxed') or (tip_type == 'r')):
        try:
            import ppafm.io as io
        except ImportError:
            raise ImportError("ppafm is required for relaxed tip scans.")
        
        Q = config['input']['Q']
        K = config['input']['K']
        data_format = config['input']['data_format']

        path_pos = os.path.join(
            config['input']['path'],
            f"Q{Q:1.2f}K{K:1.2f}/"
        )
        tip_r, lvec, nDim, atomic_info_or_head = io.load_vec_field(
            path_pos+'PPpos',
            data_format=data_format
        )

        extent = (
            lvec[0,0], lvec[0,0]+lvec[1,0],
            lvec[0,1], lvec[0,1]+lvec[2,1]
        )

        dx=lvec[1,0]/(nDim[2]-1)
        dy=lvec[2,1]/(nDim[1]-1)
        dz=lvec[3,2]/(nDim[0]-1)
        tip_r0 = RS.mkSpaceGrid(
            lvec[0,0], lvec[0,0]+lvec[1,0], dx,
            lvec[0,1], lvec[0,1]+lvec[2,1], dy,
            lvec[0,2], lvec[0,2]+lvec[3,2], dz
        )

    else:
        sw = np.asarray(config['scan']['scan_window'])
        scan_dim = np.asarray(config['scan']['scan_dim'])
        extent = (
            sw[0,0], sw[1,0],
            sw[0,1], sw[1,1],
        )
        dx = (sw[1,0]-sw[0,0])/scan_dim[0]
        dy = (sw[1,1]-sw[0,1])/scan_dim[1]
        dz = (sw[1,2]-sw[0,2])/scan_dim[2]
        tip_r = RS.mkSpaceGrid(
            sw[0,0], sw[1,0], dx,
            sw[0,1], sw[1,1], dy,
            sw[0,2], sw[1,2], dz
        )
        lvec = np.asarray([
            [sw[0,0], sw[0,1], sw[0,2]],
            [sw[1,0]-sw[0,0], 0., 0.],
            [0., sw[1,1]-sw[0,1], 0.],
            [0., 0., sw[1,2]-sw[0,2]]
        ])
        tip_r0 = tip_r

        if config['output']['plot_atoms']:
            atoms, _, _ = BU.loadAtoms(
                os.path.join(
                    config['input']['path'],
                    'input_plot.xyz'
                )
            )
            atomic_info_or_head = (atoms, lvec)
        else:
            atomic_info_or_head = (np.zeros((4,1)), np.zeros((4,3)))

    return tip_r, tip_r0, lvec, extent, atomic_info_or_head

def run_stm_scan(
    config: dict,
    eigs: np.ndarray,
    coefs: np.ndarray,
    tip_positions: np.ndarray,
    atoms: np.ndarray,
) -> Tuple:
    """
    Helper function for running the STM calculations for a given configuration.

    Args:
        config (dict): The configuration dictionary.
        eigs (np.ndarray): Eigenenergies.
        coefs (np.ndarray): KS coefficients.
        tip_positions (np.ndarray): Tip positions.
        atoms (np.ndarray): Atom positions, types, charges.

    Returns:
        result (Tuple): (STM, didv) - The result of the STM calculations.
    """

    V = config['scan']['V']
    V_max = config['scan']['V_max']
    dV = config['scan']['dV']

    scan_type = config['scan']['scan_type']
    didv_args = {
        'WF': config['advanced']['work_function'],
        'eta': config['scan']['eta'],
        'eig': eigs,
        'R': tip_positions,
        'Rat': atoms,
        'coes': coefs,
        'orbs': config['scan']['sample_orbs'],
    }
    tip_coefficients = get_tip_coefficients(config['scan']['tip_orb'])
    
    if scan_type in ['didv', 'dIdV', 'didv-single']:
        didv = PS.dIdV(
            V,
            **didv_args,
            **tip_coefficients,
        )
        # Add one dimension to the didv array
        didv = didv[None, :]
        current = None

    elif scan_type in ['states', 'STATES']:
        states = np.sort(eigs)
        mask = (states >= V) & (states <= V_max)
        states = states[mask]

        didvs = []
        for state in states:
            didv = PS.dIdV(
                V=state,
                **didv_args,
                **tip_coefficients,
            )
            didvs.append(didv)
        didv = np.stack(didvs, axis=0)
        current = None

    elif scan_type in ['STM', 'STM-single']:
        nV = abs(V/dV)+1
        current = PS.STM(
            V,
            nV,
            WF_decay=config['advanced']['work_function_decay'],
            **didv_args,
            **tip_coefficients,
        )
        # Add one dimension to the current array
        current = current[None, :]
        didv = None
    elif scan_type in ['v-scan', 'V-scan', 'Voltage-scan']:
        current, didv = PS.MSTM(
            V,
            V_max,
            dV=dV,
            WF_decay=config['advanced']['work_function_decay'],
            **didv_args,
            **tip_coefficients,
        )
    else:
        raise ValueError(f"Invalid scan type: {scan_type}")
    
    return current, didv
