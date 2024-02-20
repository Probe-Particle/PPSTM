import os
import sys
import time

import numpy as np
import h5py

from . import ProbeSTM as PS
from . import ReadSTM  as RS
from . import OCL      as ocl 

#root = os.path.abspath('../PPAFM/')
#ocl_path = "/u/79/kurkil1/unix/work/PPAFM_ocl"
#sys.path.append(ocl_path)

import ppafm as PP
import ppafm.ocl.AFMulator as afm

class STMulator():
    """
    Class for simulating scanning tunneling microscopy on molecules. 
    Requires OpenCL_py3 branch of ProbeParticleModel to be installed to the same folder.

    Functionalities:
        1. Runs AFM scan on a molecule and saves the relaxed tip positions
        2. Runs relaxed STM/didv scan on a molecule with given parameters.
        3. Returns scan result as an np.ndarray 

    Arguments:
        scan_type: string. STM or didv. (default: "didv") 
        bias_voltage: float or tuple. (V or (V_min, V_max, dV)) 
            Bias voltage to use while scanning. Tuple format required for STM scan (defeult: 0.0)
        eta: float. broadening factor of the states. (default: 0.1)
        work_function: float. (eV) (default: 5.0)
        work_function_decay: float. rate of change for work function w.r.t. bias voltage (default: 0.5)
        tip_orb: string. the tip orbital [s, pxy, spxy, 5spxy, 10spxy, CO, pz, dz2, dxzyz] (default: s)
        sample_orb : string. orbitals of the sample 'sp' (light atoms only, faster) or 'spd' (all atoms) (default: 'sp')
        scan_window: tuple ((x1, y1, z1), (x2, y2, z2)). start and end positions for scanning in each direction. 
        scan_dim: tuple (nx, ny, nz). number of tip steps in each direction.
        tip_type: string. [fixed, relaxed]. use relaxed or fixed tip. (default: 'relaxed')
        lvec: np.ndarray of shape (4, 3) or None. Unit cell boundaries for force field. First (row) vector
            specifies the origin, and the remaining three vectors specify the edge vectors of the unit cell.
            If None, will be calculated automatically from scan_window and tipR0, leaving some additional space
            on each side.
        pixPerAngstrome: int. Number of pixels (voxels) per angstrom in force field grid.
        iZPPs: int. Element of probe particle.
        qzs: Array of lenght 4. Position tip charges along z-axis relative to probe-particle center in angstroms.
        qs: Array of lenght 4. Values of tip charges in units of e/A^dim. Some can be zero.
        tip_r0: np.ndarray([x0, y0, z0]). Tip equilibrium position 
        tip_stiffness: Array of lenght 4. [N/m] (x,y,z,R). stiffness of harmonic springs anchoring probe-particle to the metalic tip-apex 
        relax_params: Array of lenght 4. (dt,damp, .., ..). parameters of relaxation, in the moment just dt and damp are used
        return_afm: boolean. Return afm scan of molecule also (default: False)
        timings: boolean. Print timing information. (default: False)
    """

    def __init__(self,
        scan_type           = "didv",
        bias_voltage        = 0.0,
        eta                 = 0.1,
        work_function       = 5.0,
        work_function_decay = 0.5,
        tip_orb             = 's',
        sample_orb          = 'sp',
        scan_window         = ((2.0, 2.0, 9.0), (18.0, 18.0, 10.0)),
        scan_dim            = (128, 128, 20),
        tip_type            = "relaxed",
        lvec                = None,
        pix_per_angstrom    = 10,
        iZPP                = 8,
        qzs                 = [0.1,  0, -0.1, 0],
        qs                  = [-10, 20,  -10, 0],
        tip_r0              = np.array([0.0, 0.0, 3.0]),
        tip_stiffness       = [0.25, 0.25, 0.0,     30.0],
        relax_params        = [0.5,  0.1,  0.1*0.2, 0.1*5.0],
        return_afm          = False,
        timings             = False
    ):
        assert scan_type in ['didv', 'STM'], "Scan type must be either didv or STM"
        assert tip_type in ['fixed', 'relaxed'], "Tip type must be either fixed or relaxed"
        if return_afm:
            assert tip_type == 'relaxed', "AFM is not performed for fixed tip STM. set return_afm=False"
        if scan_type == 'STM':
            assert isinstance(bias_voltage, tuple), "For STM scan, bias voltages must be in format (V_min, V_max, dV)."
        if scan_type == 'didv':
            assert isinstance(bias_voltage, float), "For didv scan, bias voltages must be a float."
        # TODO: more testing for input parameters

        self.scan_type = scan_type
        self.bias_voltage = bias_voltage
        self.scan_window = scan_window
        self.scan_dim = scan_dim
        self.eta = eta
        self.work_function = work_function
        self.work_function_decay = work_function_decay
        self.tip_orb = tip_orb
        self.sample_orb = sample_orb

        self.lvec = lvec
        self.pix_per_angstrom = pix_per_angstrom
        
        self.tip_type = tip_type
        self.iZPP = iZPP
        self.qzs = qzs
        self.qs = qs
        self.tip_r0 = tip_r0
        self.tip_stiffness = tip_stiffness
        self.relax_params = relax_params
        self.tip_xyz = None
        self.return_afm = return_afm
        self._get_tip_orb_coefs()

        if self.tip_type == 'fixed':   # prepare grid containing fixed tip positions, relaxed grid is calculated again for each molecule
            self._init_fixed_tip()

        if self.tip_type == 'relaxed':
            self._init_afmulator()
            #self.afmulator = AFMulator(scan_dim=self.scan_dim, scan_window=self.scan_window)

        self.timings = timings

    def eval(self, xyz, Zs, qs, eigs, coefs, REAs=None):
        """
        Prepare and evaluate STM and AFM scan.
        Arguments:
            xyz: np.ndarray(N,3). Carteesian coordinates of each atom in the system
            Zs: np.ndarray(N). Atomic numbers of each atom
            qs: np.ndarray(N). Point charges of each atom
            eigs: np.ndarray. Eigenenergies of the system
            coefs: np.ndarray. KS coefficients of the system
            REAs: np.ndarray of shape (num_atoms, 4). Lennard Jones interaction parameters. Calculated automatically if None.
        Returns:
            if return_afm==True:
                tuple of np.ndarrays. STM/didv and AFM images
            else:
                np.ndarray. STM/didv image
        """

        if REAs is None:
            self.REAs = PP.getAtomsREA(self.iZPP, Zs, self.afmulator.typeParams, alphaFac=-1.0)
        else:
            self.REAs = REAs

        if self.tip_type == 'relaxed':
            if self.timings: afm_start=time.time()
            x_afm = self.afmulator(xyz, Zs, qs, REAs=self.REAs)
            if self.timings: print(f"AFM time: {time.time()-afm_start:.2f} s")
            if self.timings: down_start=time.time()
            paths = self.afmulator.scanner.downloadPaths()
            self.tip_xyz = paths.astype(np.float64)#.transpose(1, 0, 2, 3)
            if self.timings: print(f"Tip download time: {time.time()-down_start:.2f} s")


        # if scan at HOMO -> shift HOMO energy to zero i.e. substract largest negative energy
        #homo = np.where(eigs<0, eigs, -np.inf).max()
        #lumo = np.where(eigs>0, eigs, np.inf).min()
        #print(homo, lumo)
        #eigs = eigs - lumo 

        i_min = (eigs < -2.5).sum()
        i_max = (eigs < 2.5).sum()

        eigs = eigs[i_min:i_max]
        coefs = coefs[i_min:i_max, :]

        if self.timings: stm_start=time.time()
        x = self._evalSTM(xyz, eigs, coefs)
        if self.timings: print(f"STM time: {time.time()-stm_start:.2f} s")
      
        if self.return_afm:
            out = x[:,:,::2], x_afm[:,:,:-1]
        else:
            out = x
        return out

    def _evalSTM(self, xyz, eigs, coefs):
        """
        Run STM scan for given molecule.
            Arguments:
            xyz: np.ndarray(N,3). Carteesian coordinates of each atom in the system
            eigs: np.ndarray. Eigenenergies of the system
            coefs: np.ndarray. KS coefficients of the system
        """
        if self.scan_type == "STM":
            kargs = (*self.bias_voltage, 
                     self.work_function, 
                     self.eta, 
                     eigs, 
                     self.tip_xyz,
                     xyz, 
                     coefs, 
                     self.sample_orb, 
                     *self.tip_orb_coefs)
            x_stm, x_didv = PS.MSTM(*kargs)
            x = x_stm, x_didv

        if self.scan_type == "didv":
            kargs = (self.bias_voltage, 
                     self.work_function, 
                     self.eta, 
                     eigs, 
                     self.tip_xyz,
                     xyz, 
                     coefs, 
                     self.sample_orb, 
                     *self.tip_orb_coefs)
            x = PS.dIdV(*kargs)
        return x 

    def ocl_stm(self, xyz, qs, coefs, eigs):
        atoms = ocl.xyzq2float4(xyzs=xyz, qs=qs)
        CAOs = ocl.CAOsp2f4(coefs, xyz.shape[0])
        spectral = ocl.getSpectral(eigs, Wf=self.work_function, w=self.bias_voltage)
        kargs = ocl.initArgs(atoms, CAOs, spectral, self.tip_xyz)
        x = ocl.run(kargs, self.tip_xyz.shape[:3])
        return x
        
    def _get_tip_orb_coefs(self):
        """Get tip orbital coefficients."""
        if (self.tip_orb == 's'):
            tc = (1.,0.,0.,0.,0.,0.,0.)         # [s, px, py, pz, dz2, dxz, dyz ] 
        elif (self.tip_orb == 'pxy'):
            tc = (0.,0.5,0.5,0.,0.,0.,0.)       # [s, px, py, pz, dz2, dxz, dyz ] 
        elif (self.tip_orb == 'spxy'):
            tc = (0.5,0.25,0.25,0.,0.,0.,0.)    # [s, px, py, pz, dz2, dxz, dyz ] 
        elif (self.tip_orb == '5spxy'):
            tc = (0.05,0.475,0.475,0.,0.,0.,0.) # [s, px, py, pz, dz2, dxz, dyz ] 
        elif (self.tip_orb == '10spxy'):
            tc = (0.10,0.45,0.45,0.,0.,0.,0.)   # [s, px, py, pz, dz2, dxz, dyz ] 
        elif (self.tip_orb == 'CO'):
            tc = (0.15,0.5,0.5,0.,0.,0.,0.)     # [s, px, py, pz, dz2, dxz, dyz ] 
        elif (self.tip_orb == 'pz'):
            tc = (0.,0.,0.,1.,0.,0.,0.)         # [s, px, py, pz, dz2, dxz, dyz ] 
        elif (self.tip_orb == 'dz2'):
            tc = (0.,0.,0.,0.,1.,0.,0.)         # [s, px, py, pz, dz2, dxz, dyz ] 
        elif (self.tip_orb == 'dxzyz'):
            tc = (0.,0.,0.,0.,0.,0.5,0.5)       # [s, px, py, pz, dz2, dxz, dyz ] 
        self.tip_orb_coefs = tc

    def _init_fixed_tip(self):
        """Initialize grid for fixed tip scan."""
        dx = (self.scan_window[1][0]-self.scan_window[0][0])/self.scan_dim[0]
        dy = (self.scan_window[1][1]-self.scan_window[0][1])/self.scan_dim[1]
        dz = (self.scan_window[1][2]-self.scan_window[0][2])/self.scan_dim[2]
        
        tip_r  = RS.mkSpaceGrid(self.scan_window[0][0],self.scan_window[1][0],dx,
                                self.scan_window[0][1],self.scan_window[1][1],dy,
                                self.scan_window[0][2],self.scan_window[1][2],dz)

        self.tip_xyz = tip_r

    def _init_afmulator(self):
        """Initializes the AFM simulator using same parameters as the STM simulator."""
        self.afmulator = afm.AFMulator(scan_dim=self.scan_dim, 
                                   scan_window=self.scan_window,
                                   lvec=self.lvec,
                                   pixPerAngstrome=self.pix_per_angstrom,
                                   iZPP=self.iZPP,
                                   Qs=self.qs,
                                   QZs=self.qzs,
                                   tipR0=self.tip_r0,
                                   tipStiffness=self.tip_stiffness)

    def __call__(self, xyz, Zs, qs, eigs, coefs, REAs=None):
        """
        Make object callable.
        Arguments:
            xyz: np.ndarray(N,3). Carteesian coordinates of each atom in the system
            Zs: np.ndarray(N). Atomic numbers of each atom
            qs: np.ndarray(N). Point charges of each atom
            eigs: np.ndarray. Eigenenergies of the system
            coefs: np.ndarray. KS coefficients of the system
            REAs: np.ndarray of shape (num_atoms, 4). Lennard Jones interaction parameters. Calculated automatically if None.
        """
        return self.eval(xyz, Zs, qs, eigs, coefs, REAs)
