##########################################################################################################################
#                                                                                                                        #
#                                  What follows are options of the PP-STM (dIdV) code                                    #
#                                                                                                                        #
##########################################################################################################################
#
# Note : This type of simulations works for solid slabs or molecules on slabs (substrate) ; for freestanding molecule it can give you nonsences
#
# ***** System information: *****

# Converting 1D and 2D string arrays to arrays of floats
import re

def conv1Darray(array1D):
    array1D = re.sub('\[|\]', '', array1D)
    newArray = array1D.split(',')
    return list(map(float, newArray))

def conv2Darray(array2d):
    newArray = array2d.split(']')
    numarray = []
    for n in newArray:
        if len(n) == 0:
            break
        n = re.sub('\[', '', n)
        n = n.split(',')
        for i in n:
            if len(i) == 0:
                n.remove(i)
        numarray.append(list(map(float, n)))
    return numarray
        
def importData(myDict, paths):

    files_path = paths['inputPath']            # where are files fron DFT code ; rather do not use this #
    
    try:
        # None ; [[ax,ay,0],[bx,by,0]],[0,0,cz]] or [[ax,ay],[bx,by]] ; 'input.lvs' -- files with specified cell ; in FHI-AIMS & GPAW allready specified with geometry #
        lvs = myDict['lvs']
        if lvs == 'None' or len(lvs) == 0: lvs = None
        elif lvs[0] == '[':
            lvs = conv2Darray(lvs)
        else:
            lvs = files_path + lvs
    except:
        print('lvs input not correct')

    # E.G. 'input.xyz' , 'input.bas' , 'geometry.in'; None for GPAW #
    geometry_file = paths['geometry_file']

    # 'fireball'='Fireball'='FIREBALL' ; 'aims'='AIMS'='FHI-AIMS' ; 'cp2k'='CP2K' ; 'gpaw'='GPAW' #
    dft_code = myDict['dft_code']

    pbc = (int(myDict['pbc'][0]), int(myDict['pbc'][1]))        # (0,0) = None = False -- only original geometry ; (0.5,0.5) -- 2x2 cell ; (1,1) -- 3x3 cell (around original) ; (2,2) -- 5x5 cell (around original) ... #
    cp2k_name = paths['cp2kName']  # Name used in CP2K calculations or GPAW calc #
    if cp2k_name == 'none': cp2k_name = None
    
    cut_atoms = int(myDict['cut_atoms'])         # None = -1 -- All atoms of the sample contributes to tunelling ; 1 -- only 1st atom of the sample contributes to the tunelling ; 57 -- first 57 atoms of the sample contributes to the tunelling ; ... #
    
    try:
        lower_atoms = myDict['lower_atoms']             # [] = None -- No atoms has lowered hopping ; be aware python numbering occurs here: [0] - means lowering of the 1st atom; [0,1,2,3] -- lowering of 1st 4 atoms ... #
        if lower_atoms == 'None' or len(lower_atoms) == 0: lower_atoms = []
        else:
            lower_atoms = list(map(int,conv1Darray(lower_atoms)))
    except:
        print('lower Atoms input not corrct')
    
    try:
        lower_coefs = myDict['lower_coefs']             # [] = None -- No lowering of the hoppings  ; [0.5] -- lowering of the 1st atom hopping to 0.5                           ; [0.5,0.5,0.5,0.5] -- lowering of 1st 4 atoms to 0.5 ... #
        if lower_coefs == 'None' or len(lower_coefs) == 0: lower_coefs = []
        else:
            lower_coefs = conv1Darray(lower_coefs)
    except:
        print('lower_coefs input not correct')

    # None=0.0 -- no change to the Fermi Level ; -0.1 -- shifts the Fermi Level by 0.1 eV lower ... #
    fermi = None
    # cut out all orbitals lower than  -2.5 eV bellow Rermi (should be: cut_min <= Vmin-2*eta) . taken to the Fermi Level #
    cut_min = -2.5
    # cut out all orbitals higher than -2.5 eV above  Fermi (should be: cut_max >= Vmax+2*eta) . taken to the Fermi Level #
    cut_max = +2.5
    sample_orbs = myDict['sample_orbs']
    spin = myDict['spin']
    if spin == 'None': spin = None
    print("Reading electronic & geometry structure files")
    import numpy as np

    if ((dft_code == 'fireball') or (dft_code == 'Fireball') or (dft_code == 'FIREBALL') or (dft_code == 'cp2k') or (dft_code == 'CP2K')):
        if isinstance(lvs, (list, tuple, np.ndarray)):
            cell = np.array([[lvs[0][0], lvs[0][1], 0.0], [lvs[1][0], lvs[1][1], 0.0], [
                                    0.0, 0.0, 99.9]]) if len(lvs) == 2 else lvs
        elif isinstance(lvs, (str)):
            cell = np.loadtxt(lvs)
        elif ((pbc == (0, 0)) or (pbc == (0., 0.))):
                    cell = np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]])
        else:
            print("PBC required, but lattice vector not specified. What can I do with that? I rather go to eat something.")
            return None
    
    from . import ReadSTM as RS
            
    eigEn = coefs = Ratin = None
    if ((dft_code == 'fireball') or (dft_code == 'Fireball') or (dft_code == 'FIREBALL')):
            eigEn, coefs, Ratin = RS.read_FIREBALL_all(name=files_path + 'phik_0001_', geom=files_path+geometry_file, lvs=cell, fermi=fermi,
                                                orbs=sample_orbs, pbc=pbc, cut_min=cut_min, cut_max=cut_max, cut_at=cut_atoms, lower_atoms=lower_atoms, lower_coefs=lower_coefs)

    elif ((dft_code == 'gpaw') or (dft_code == 'GPAW')):
            eigEn, coefs, Ratin = RS.read_GPAW_all(name=files_path + cp2k_name + '.gpw', fermi=fermi, orbs=sample_orbs,
                                            pbc=pbc, cut_min=cut_min, cut_max=cut_max, cut_at=cut_atoms, lower_atoms=lower_atoms, lower_coefs=lower_coefs);

    elif ((dft_code == 'aims') or (dft_code == 'AIMS') or (dft_code == 'FHI-AIMS')):
        if ((spin == None) or (spin == False)):
            name = 'KS_eigenvectors.band_1.kpt_1.out'
        elif ((spin == 'up') or (spin == 'alpha') or (spin == 'both')):
                name = 'KS_eigenvectors_up.band_1.kpt_1.out'
        elif ((spin == 'down') or (spin == 'beta') or (spin == 'dn')):
                    name = 'KS_eigenvectors_dn.band_1.kpt_1.out'
        else:
            print("unknown spin, I'm going to sleep. Good Night")
            return None

        eigEn, coefs, Ratin = RS.read_AIMS_all(name=files_path + name, geom=files_path + geometry_file, fermi=fermi, orbs=sample_orbs,
                                            pbc=pbc, cut_min=cut_min, cut_max=cut_max, cut_at=cut_atoms, lower_atoms=lower_atoms, lower_coefs=lower_coefs);
        if (spin == 'both'):
            eigEn1 = eigEn.copy(); coefs1 = coefs.copy(); del eigEn, coefs;
            name = 'KS_eigenvectors_dn.band_1.kpt_1.out'
            eigEn2, coefs2, Ratin = RS.read_AIMS_all(name=files_path + name, geom=files_path + geometry_file, fermi=fermi, orbs=sample_orbs,
                                                    pbc=pbc, cut_min=cut_min, cut_max=cut_max, cut_at=cut_atoms, lower_atoms=lower_atoms, lower_coefs=lower_coefs);
            eigEn = np.concatenate((eigEn1, eigEn2), axis=0)
            coefs = np.concatenate((coefs1, coefs2), axis=0)

    elif ((dft_code == 'cp2k') or (dft_code == 'CP2K')):
        if ((spin == None) or (spin == False)):
            eigEn, coefs, Ratin = RS.read_CP2K_all(name=files_path + cp2k_name, lvs=cell, fermi=fermi, orbs=sample_orbs,
                                                pbc=pbc, cut_min=cut_min, cut_max=cut_max, cut_at=cut_atoms, lower_atoms=lower_atoms, lower_coefs=lower_coefs);
        elif ((spin == 'up') or (spin == 'alpha')):
            eigEn, coefs, Ratin = RS.read_CP2K_all(name=files_path + cp2k_name, lvs=cell, fermi=fermi, orbs=sample_orbs, pbc=pbc,
                                                cut_min=cut_min, cut_max=cut_max, cut_at=cut_atoms, lower_atoms=lower_atoms, lower_coefs=lower_coefs, spin='alpha');
        elif (spin == 'both'):
            eigEn1, coefs1, Ratin = RS.read_CP2K_all(name=files_path + cp2k_name, lvs=cell, fermi=fermi, orbs=sample_orbs, pbc=pbc,
                                                    cut_min=cut_min, cut_max=cut_max, cut_at=cut_atoms, lower_atoms=lower_atoms, lower_coefs=lower_coefs, spin='alpha');
            eigEn2, coefs2, Ratin = RS.read_CP2K_all(name=files_path + cp2k_name, lvs=cell, fermi=fermi, orbs=sample_orbs, pbc=pbc,
                                                    cut_min=cut_min, cut_max=cut_max, cut_at=cut_atoms, lower_atoms=lower_atoms, lower_coefs=lower_coefs, spin='beta');
            eigEn = np.concatenate((eigEn1, eigEn2), axis=0)
            coefs = np.concatenate((coefs1, coefs2), axis=0)
        elif ((spin == 'down') or (spin == 'beta') or (spin == 'dn')):
            eigEn, coefs, Ratin = RS.read_CP2K_all(name=files_path + cp2k_name, lvs=cell, fermi=fermi, orbs=sample_orbs, pbc=pbc,
                                                cut_min=cut_min, cut_max=cut_max, cut_at=cut_atoms, lower_atoms=lower_atoms, lower_coefs=lower_coefs, spin='beta');
        else:
            print("unknown spin, I'm going to sleep. Good Night")
            return None
    
    return {'eigEn': eigEn, 'coefs': coefs, 'Ratin': Ratin}

    
def conv2float(string):
    try:
        return float(string)
    except ValueError:
        return string

def newPPSTM_simple(myDict, paths, importData):
    
    eigEn = importData['eigEn']
    coefs = importData['coefs']
    Ratin = importData['Ratin']

    tip_type = myDict['tip_type']
    #
    ncpu = myDict['OMP_NUM_THREADS']               # number of cpu cores for OMP paralelization: ncpu = 1 -- serial compilation & calculations; iff ncpu > 1, then OMP paralel recompilation is used and C++ calculations are running on more cores #
    #
    # ***** Main informations ******
    #
    scan_type = myDict['scan_type']     # 'didv'='dIdV''='didv-single' -- only dIdV for one voltage = V ; 'v-scan'='V-scan'='Voltage-scan' -- both STM & dIdV scan - V .. Vmax; 'STM'='STM-single' -- STM for one Voltage = V, use V-scan rather #
    # 'relaxed'    # 'fixed'='f' -- for stiff/metal tip apexes ; 'relaxed'='r' -- for flexible tip apexes (precalculated by PP-AFM) . For this option you have to have "installed" PPAFM in your PPSTM directory #
    
    V = myDict['V']         # !!!! V = Vmin for SCAN !!!! #
    # V = V_min >= -2.0 V ; V_max <= 2.0 V (othervise changes in the later code needed) #
    V_max = myDict['Vmax']
    dV = myDict['dV']      # voltage step , dV <= 0.1 V #
    # Lorentzian width of states in energy scale: typically 0.1; can be in range of 0.3-0.05 eV in some cases (low amount of layers ...) even up to 1.0 eV #
    eta = myDict['etaValue']
    # 0.0 <= WF_decay <= 1.0 ; How fast WorkFunction tunnelling barrier is changing with Voltage : (WF = WF_0 + V*WF_decay) -- 0.0 no change ; 1.0 - the same change as voltage #
    WF_decay = myDict['wf_decay']
    # orbitals of the sample 'sp' (light atoms only, faster) or 'spd' (all atoms) #
    sample_orbs = myDict['sample_orbs']

    # None=False ; for FHI-AIMS & CP2K: None -- spin-unpolarized/spin-restricted calc. ;  'both' , 'up'='alpha' or 'down" (last 3 are for spin-polarizes or spin-unrestricted calculations) #
    spin = myDict['spin']
    if spin == 'None': spin = None
    #
    # ***** Informations for x,y,z tip_type = 'fixed' ******
    #
    x = myDict['x']     # [xmin, xmax, dx] #
    y = myDict['y']     # [ymin, ymax, dy] #
    # !!!! z-starts from zero - normally zmin >= 3-4 Ang above heighest atoms !!!! [zmin, zmax, dz] ; for single height scan use : [z, z, 1.0] #
    z = myDict['z']
    #
    # ***** Informations for PP positions, tip_type = 'relaxed' ******
    #
    # charge (PP-AFM) ; Ocharge PP-AFM complex_tip autumn 2018) ; [e] (monopole), [e*A] (dipole), [e*A^2] (quadrupole) #
    Q = myDict['qValue']
    # x stiffness (PP-AFM master autumn 2018); klat (PP-AFM dev/OpenCl autumn 2018); Oklat (PP-AFM complex_tip autumn 2018) ; [N/m] #
    K = myDict['kValue']
    # 'xsf'='XSF' ; 'npy'='NPY' ; -- format in which PPpos are stored from PP-AFM run #
    data_format = myDict['data_format']
    #
    # *****Output options ******
    #
    # True / False -- plot "png" images (2D constant height) #
    # PNG = True
    # True / False -- write ".xyz" WSxM files (2D constant height) #
    # WSxM = False
    # XSF = False                 # True / False -- write ".xsf" files with 3D stucks of data . For this option you have to have "installed" PPAFM in your PPSTM directory #
    # NPY = False                 # True / False -- write ".npy" numpy binary files with 3D stucks of data . For this option you have to have "installed" PPAFM in your PPSTM directory #
    # True / False -- plot geometry (position of atoms into the PNG images and into the XSF files). You have to have your geometry, which you want to plot in input_plot.xyz. This doesn't change the name of the output files #
    # plot_atoms = True
    WorkFunction = 5.0          # 5.0 eV is standart #
    #
    #
    ##########################################################################################################################
    #                                                                                                                        #
    #                                 DO NOT TOUCH LATER CODE (unless you know what you are doing)                           #
    #                                                                                                                        #
    ##########################################################################################################################

    print("Importing libraries")
    
    import os
    from . import ReadSTM as RS
    from . import ProbeSTM as PS
    import numpy as np
    import sys
    sys.path.append('../') 

    if (ncpu > 1):
        os.environ['OMP_NUM_THREADS'] = str(ncpu)
        print('OMP_NUM_THREADS:', os.environ['OMP_NUM_THREADS'])

    if (tip_type == 'relaxed') or (tip_type == 'r'):
        print("For XSF or NPY outputs or tip_type = relaxed you have to have installed PPAFM in your PPSTM directory ")
        from . import GridUtils as GU
    
    print("Libraries imported")
    # --- Initial check --- #

    # assert(PNG or WSxM or XSF or NPY), "No output set to be True; I'm not going to do anything if there is no output. I'm too lazy like a Gartfield. "

    # --- specification of tip orbitals --- #
    # 's' ; 'pxy' -- px & py ; 'spxy' -- 50% s & 50% pxy ; '5spxy' -- 5% s & 95% pxy ; '10spxy' -- 10% s & 90% pxy ; 'CO' -- 13% s & 87% pxy (PRL 119, 166001 (2017)) ; 'pz' ; For sample_orbs = 'sp' , possible 'dz2' and 'dxzyz' -- dxz & dyz #
    s = myDict['tipOrbS']
    px = py = 0.5 * myDict['tipOrbPxy']
    tc = [s, px, py, 0., 0., 0., 0.]

    # if (tip_orb == 's'):
    #  tc = [s, px, py, 0., 0., 0., 0.]  # [s, px, py, pz, dz2, dxz, dyz ]
    # elif (tip_orb == 'pxy'):
    #  tc = [s, px, py, 0., 0., 0., 0.]  # [s, px, py, pz, dz2, dxz, dyz ]
    # elif (tip_orb == 'spxy'):
    #  tc = [s, px, py, 0., 0., 0., 0.]  # [s, px, py, pz, dz2, dxz, dyz ]
    # elif (tip_orb == '5spxy'):
        # [s, px, py, pz, dz2, dxz, dyz ]
    #  tc = [s, px, py, 0., 0., 0., 0.]
    # elif (tip_orb == '10spxy'):
    #  tc = [s, px, py, 0., 0., 0., 0.]  # [s, px, py, pz, dz2, dxz, dyz ]
    # elif (tip_orb == 'CO'):
    #  tc = [s, px, py, 0., 0., 0., 0.]  # [s, px, py, pz, dz2, dxz, dyz ]
    # elif (tip_orb == 'pz'):
    #  tc = [s, px, py, 1., 0., 0., 0.]  # [s, px, py, pz, dz2, dxz, dyz ]
    # elif (tip_orb == 'dz2'):
    #  tc = [s, px, py, 0., 1., 0., 0.]  # [s, px, py, pz, dz2, dxz, dyz ]
    # elif (tip_orb == 'dxzyz'):
    #  tc = [s, px, py, 0., 0., 0.5, 0.5]  # [s, px, py, pz, dz2, dxz, dyz ]
    # else:
    #  print("Don't know what kind of tip you mean. I rather going to exit."); return None

    # print "DEBUG: tc ", tc , " [s, px, py, pz, dz2, dxz, dyz ] "

    # --- the grid on which the STM signal is calculated --- #

    if tip_type == 'relaxed':
        try:
            print("Importing positions of PP from the PP-AFM calculations. Path for the data:")
            path_pos = "Q%1.2fK%1.2f/" % (Q, K)
            print(path_pos)
            tip_r, lvec, nDim = GU.load_vec_field(
                os.path.join(paths['inputPath'], path_pos+'PPpos'), data_format=data_format)
            extent = (lvec[0, 0], lvec[0, 0]+lvec[1, 0],
                    lvec[0, 1], lvec[0, 1]+lvec[2, 1])
            # print "DEBUG: extent", extent
            print("PP postions imported")
            dx = lvec[1, 0]/(nDim[2]-1); dy = lvec[2, 1] / \
                            (nDim[1]-1); dz = lvec[3, 2]/(nDim[0]-1);
            tip_r0 = RS.mkSpaceGrid(lvec[0, 0], lvec[0, 0]+lvec[1, 0], dx, lvec[0, 1],
                                    lvec[0, 1]+lvec[2, 1], dy, lvec[0, 2], lvec[0, 2]+lvec[3, 2], dz)
            # print "DEBUG: dx, dy, dz", dx, dy, dz
            # print "DEBUG: tip_r.shape, tip_r0.shape", tip_r.shape, tip_r0.shape
        except:
            print('Relaxed scan not possible. Firstly you neeed to install PPAFM code and run < ./run_test to create pre-calculated positions.')
            return None
    else:
        print("Priparing the scan grid for fixed scan")
        extent = (x[0], x[1], y[0], y[1])
        tip_r = RS.mkSpaceGrid(x[0], x[1], x[2], y[0],
                            y[1], y[2], z[0], z[1], z[2])
        lvec = np.array([[x[0], y[0], z[0]], [x[1]-x[0], 0., 0.],
                        [0., y[1]-y[0], 0.], [0., 0., z[1]-z[0]]])
        # print "DEBUG: extent", extent
        # print "DEBUG: lvec", lvec
        tip_r0 = tip_r
        print("scan grids prepared")
    
    # --- the Main calculations --- #
    # 'didv'='dIdV''='didv-single' -- only dIdV for one voltage = V ; 'v-scan'='V-scan'='Voltage-scan' -- both STM & dIdV scan - V .. Vmax; 'STM'='STM-single' -- STM for one Voltage = V, use V-scan rather #

    states_b = False

    if ((scan_type == 'didv') or (scan_type == 'dIdV') or (scan_type == 'didv-single')):
        didv = np.array([PS.dIdV(V,    WorkFunction, eta, eigEn, tip_r, Ratin, coefs, orbs=sample_orbs,
                        s=tc[0], px=tc[1], py=tc[2], pz=tc[3], dz2=tc[4], dxz=tc[5], dyz=tc[6])])
        WF_decay = 0.0;
        # print "DEBUG: didv.shape ", didv.shape
    elif ((scan_type == 'states') or (scan_type == 'STATES')):
        states = np.sort(eigEn); mask = states >= V; states = states[mask]; del mask;
        mask = states <= V_max; states = states[mask]; del mask;
        fst = True
        print("DEBUG: states:", states)
        for isi in states:
            if fst:
                didv = np.array([PS.dIdV(isi, WorkFunction, eta, eigEn, tip_r, Ratin, coefs, orbs=sample_orbs,
                                s=tc[0], px=tc[1], py=tc[2], pz=tc[3], dz2=tc[4], dxz=tc[5], dyz=tc[6])])
                fst = False
            else:
                didv = np.append(didv, [PS.dIdV(isi,    WorkFunction, eta, eigEn, tip_r, Ratin, coefs, orbs=sample_orbs,
                                s=tc[0], px=tc[1], py=tc[2], pz=tc[3], dz2=tc[4], dxz=tc[5], dyz=tc[6])], axis=0)
        # print "DEBUG: didv.shape ", didv.shape
        states_b = True; WF_decay= 0.0;
    elif ( (scan_type == 'STM') or (scan_type == 'STM-single') ):
        nV = abs(V/dV)+1
        # print "DEBUG: V, nV:", V, nV
        current = np.array([   PS.STM( V, nV, WorkFunction, eta, eigEn, tip_r, Ratin, coefs, orbs=sample_orbs, s=tc[0], px =tc[1], py=tc[2], pz=tc[3], dz2=tc[4], dxz=tc[5], dyz=tc[6], WF_decay=WF_decay) ])
        
        # print "DEBUG: current.shape ", current.shape
    else:
        current, didv = PS.MSTM( V, V_max, dV, WorkFunction, eta, eigEn, tip_r, Ratin, coefs, orbs=sample_orbs, s=tc[0], px =tc[1], py=tc[2], pz=tc[3], dz2=tc[4], dxz=tc[5], dyz=tc[6], WF_decay=WF_decay)
        # print "DEBUG: didv.shape ", didv.shape
        # print "DEBUG: current.shape ", current.shape

    # =========== Utils for plotting atoms =========================


    # --- plotting part here, plots all calculated signals --- #

    Voltages= np.arange(V,V_max+0.001,dV) if not states_b else states # this part is important for scans over slabs at different voltages
    round_index = 2 if not states_b else 5
    print("Voltages", Voltages)
    namez = []
    for V in Voltages:
        namez.append(str(round(V,round_index)))
        
    plotData = {}

    try:
        NoV_didv = len(didv) 
        NoH_didv = len(didv[0]) 
        plotData['didv'] = didv
        plotData['NoV_didv'] = NoV_didv
        plotData['NoH_didv'] = NoH_didv
    except:
        plotData['didv'] = None

    try:
        NoV_STM = len(current) 
        NoH_STM = len(current[0]) 
        plotData['current'] = current
        plotData['NoV_STM'] = NoV_STM
        plotData['NoH_STM'] = NoH_STM
    except:
        plotData['current'] = None

    print() 
    print()
    print("Done")
    print()

    plotData.update({'namez': namez,
                'tip_type': tip_type,
                'lvec': lvec,
                'extent': extent,
                'WorkFunction': WorkFunction,
                'Voltages': Voltages,
                'WF_decay': WF_decay,
                'tip_r0': tip_r0,
                'tip_r': tip_r})

    return plotData
    # print "DEBUG: Voltages", Voltages
    # print "DEBUG: namez", namez
    # print "DEBUG: NoV", NoV
    # print "DEBUG: NoH", NoH
    
    # --- the end --- #

