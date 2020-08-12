#!/usr/bin/python
#
##########################################################################################################################
#                                                                                                                        #
#                         What follows are options for summing tip-orbitals' contributions                               #
#                                                                                                                        #
##########################################################################################################################
#
# Note : This type of simulations works for solid slabs or molecules on slabs (substrate) ; for freestanding molecule it can give you nonsences
#
# ***** System information: path (absolute or relative) to your PPSTM code *****
#
ppstm_path = './PPSTM/'
#
# ***** Main informations ******
#
scan_type     = 'V-scan'       # 'didv'='dIdV''='didv-single' -- only dIdV for one voltage = V ; 'v-scan'='V-scan'='Voltage-scan' -- both STM & dIdV scan - V .. Vmax; 'STM'='STM-single' -- STM for one Voltage = V, use V-scan rather #
tip_type      = 'relaxed'    # 'fixed'='f' -- for stiff/metal tip apexes ; 'relaxed'='r' -- for flexible tip apexes (precalculated by PP-AFM) . For this option you have to have "installed" PPAFM in your PPSTM directory #
V             = -0.5         # !!!! V = Vmin for SCAN !!!! #
V_max         = +0.0         # V = V_min >= -2.0 V ; V_max <= 2.0 V (othervise changes in the later code needed) #
dV            =  0.1         # voltage step , dV <= 0.1 V #
eta           =  0.1         # Lorentzian width of states in energy scale: typically 0.1; can be in range of 0.3-0.05 eV in some cases (low amount of layers ...) even up to 1.0 eV #
WF_decay      =  0.5         # 0.0 <= WF_decay <= 1.0 ; How fast WorkFunction tunnelling barrier is changing with Voltage : (WF = WF_0 + V*WF_decay) -- 0.0 no change ; 1.0 - the same change as voltage #
tip_orb1      = 's'          # 's' ; 'pxy' -- px & py ; 'spxy' -- 50% s & 50% pxy ; '5spxy' -- 5% s & 95% pxy ; '10spxy' -- 10% s & 90% pxy ; 'CO' -- 13% s & 87% pxy (PRL 119, 166001 (2017)) ; 'pz' ; For sample_orbs = 'sp' , possible 'dz2' and 'dxzyz' -- dxz & dyz #
tip_orb2      = 'pxy'        # 's' ; 'pxy' -- px & py ; 'spxy' -- 50% s & 50% pxy ; '5spxy' -- 5% s & 95% pxy ; '10spxy' -- 10% s & 90% pxy ; 'CO' -- 13% s & 87% pxy (PRL 119, 166001 (2017)) ; 'pz' ; For sample_orbs = 'sp' , possible 'dz2' and 'dxzyz' -- dxz & dyz #
data_format   = 'xsf'        # 'xsf'='XSF' ; 'npy'='NPY' ; -- format in which precalculated date are stored from 1st two PPSTM runs #
#
# *****Output options ******
#
tip_orb1_amount = 0.50       #
tip_orb2_amount = 0.50       # 
PNG  = True                  # True / False -- plot "png" images (2D constant height) #
WSxM = False                 # True / False -- write ".xyz" WSxM files (2D constant height) #
XSF  = False                 # True / False -- write ".xsf" files with 3D stucks of data . For this option you have to have "installed" PPAFM in your PPSTM directory #
NPY  = False                 # True / False -- write ".npy" numpy binary files with 3D stucks of data . For this option you have to have "installed" PPAFM in your PPSTM directory #
plot_atoms = True            # True / False -- plot geometry (position of atoms into the PNG images and into the XSF files). You have to have your geometry, which you want to plot in input_plot.xyz. This doesn't change the name of the output files #
#
# ***** More advanced options ******
#
WorkFunction =  5.0          # 5.0 eV is standart #
files_path   = ''            # where are files from 1st two PPSTM runs ; rather do not use this #
#
#
##########################################################################################################################
#                                                                                                                        #
#                                 DO NOT TOUCH LATER CODE (unless you know what you are doing)                           #
#                                                                                                                        #
##########################################################################################################################

print("Importing libraries")

import os
import sys
sys.path.append(ppstm_path) 

import numpy as np
import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend. ## !!! important for working on clusters !!!!
import matplotlib.pyplot as plt
print("For XSF or NPY outputs and inputs you have to have installed PPAFM in your PPSTM directory ")
import pyPPSTM.ReadSTM           as RS
import pyProbeParticle.GridUtils as GU
if (plot_atoms):
    import pyPPSTM.basUtils as Bu
    import pyPPSTM.elements as elements

print("Libraries imported")

# --- Initial check --- #

assert( PNG or WSxM or XSF or NPY ), "No output set to be True; I'm not going to do anything if there is no output. I'm too lazy like a Gartfield. "
didv_b   = False
STM_b    = False
V_scan_b = False

if ( (scan_type == 'didv') or (scan_type == 'dIdV') or (scan_type == 'didv-single')):
    didv_b = True; WF_decay = 0.0;
elif ( (scan_type == 'STM') or (scan_type == 'STM-single') ):
    STM_b = True
else:
    didv_b   = True
    STM_b    = True
    V_scan_b = True

# --- Reading geometry for plotting (if necessary)  --- #

if (plot_atoms):
    geom_plot, tmp1, tmp2 = Bu.loadAtoms('input_plot.xyz'); del tmp1, tmp2;
else:
    geom_plot = None
#print "DEBUG: geom_plot", geom_plot

# --- importing STM (dIdV) calculated signal on the grid --- #

Voltages=np.arange(V,V_max+0.001,dV) if V_scan_b else [V] # this part is important for scans over slabs at different voltages
namez = []
for V in Voltages:
    namez.append(str(round(V,2)))

#print "DEBUG: Voltages", Voltages

for i in range(len(Voltages)):
    if didv_b :
        print("Importing dIdV data for V:", namez[i])
        #print("DEBUG:WF:", WorkFunction , "Voltages[i]:", Voltages[i], "WF_decay", WF_decay  )
        name_file1 =  'didv_'+namez[i]+"_tip_"+tip_type+"-"+tip_orb1+"_WF_"+str(WorkFunction-Voltages[i]*WF_decay)+"_eta_"+str(eta) # WorkFunction gets higher with lower (higher negative) sample bias
        name_file2 =  'didv_'+namez[i]+"_tip_"+tip_type+"-"+tip_orb2+"_WF_"+str(WorkFunction-Voltages[i]*WF_decay)+"_eta_"+str(eta)
        #print ("DEBUG: name_file1", name_file1)
        tmp_dIdV1, lvec1, nDim1 = GU.load_scal_field( files_path+name_file1 ,data_format=data_format)
        didv1 =  np.array([tmp_dIdV1]) if i == 0 else np.append(didv1, np.array([tmp_dIdV1]),axis=0)
        #print "DEBUG: name_file2", name_file2
        tmp_dIdV2, lvec2, nDim2 = GU.load_scal_field( files_path+name_file2 ,data_format=data_format)
        didv2 =  np.array([tmp_dIdV2]) if i == 0 else np.append(didv2, np.array([tmp_dIdV2]),axis=0)
        assert np.array(lvec1).all() == np.array(lvec2).all(), "lvec1 != lvec2 control your input files"; assert np.array(nDim2).all() == np.array(nDim2).all(), "nDim1 != nDim2 control your input files"
        print("dIdV for V:",namez[i]," imported")
        #print "DEBUG: didv1.shape", didv1.shape
    if STM_b :
        print("Importing STM data for V:", namez[i])
        name_file1 =  'STM_'+namez[i]+"_tip_"+tip_type+"-"+tip_orb1+"_WF_"+str(WorkFunction)+"_WF_decay_"+str(round(WF_decay,1))+"_eta_"+str(eta)
        name_file2 =  'STM_'+namez[i]+"_tip_"+tip_type+"-"+tip_orb2+"_WF_"+str(WorkFunction)+"_WF_decay_"+str(round(WF_decay,1))+"_eta_"+str(eta)
        #print "DEBUG: name_file1", name_file1
        tmp_stm1, lvec1, nDim1 = GU.load_scal_field( files_path+name_file1 ,data_format=data_format)
        current1 =  np.array([tmp_stm1]) if i == 0 else np.append(current1, np.array([tmp_stm1]),axis=0)
        #print "DEBUG: name_file2", name_file2
        tmp_stm2, lvec2, nDim2 = GU.load_scal_field( files_path+name_file2 ,data_format=data_format)
        current2 =  np.array([tmp_stm2]) if i == 0 else np.append(current2, np.array([tmp_stm2]),axis=0)
        assert np.array(lvec1).all() == np.array(lvec2).all(), "lvec1 != lvec2 control your input files"; assert np.array(nDim2).all() == np.array(nDim2).all(), "nDim1 != nDim2 control your input files"
        print("STM for V:",namez[i]," imported")
        #print "DEBUG: current1.shape", current1.shape
lvec = lvec1; nDim = nDim1;
extent = (lvec[0,0],lvec[0,0]+lvec[1,0],lvec[0,1],lvec[0,1]+lvec[2,1])
#print "DEBUG: extent", extent
dx=lvec[1,0]/(nDim[2]-1); dy=lvec[2,1]/(nDim[1]-1); dz=lvec[3,2]/(nDim[0]-1);
tip_r0 = RS.mkSpaceGrid(lvec[0,0],lvec[0,0]+lvec[1,0],dx,lvec[0,1],lvec[0,1]+lvec[2,1],dy,lvec[0,2],lvec[0,2]+lvec[3,2],dz)
#print "DEBUG: tip_r0", tip_r0

# --- main part --- #

if STM_b :
    current = tip_orb1_amount * current1 + tip_orb2_amount * current2
if didv_b :
    didv    = tip_orb1_amount * didv1    + tip_orb2_amount * didv2

# =========== Utils for plotting atoms =========================

def plotAtoms( atoms, atomSize=0.1, edge=True, ec='k', color='w' ):
    plt.fig = plt.gcf()
    es = atoms[0]; xs = atoms[1]; ys = atoms[2]
    for i in range(len(xs)):
        fc = '#%02x%02x%02x' % elements.ELEMENT_DICT[es[i]][7] #; print "DEBUG: fc", fc ; ##fc = '#FFFFFF' ##
        if not edge:
            ec=fc
        circle=plt.Circle( ( xs[i], ys[i] ), atomSize, fc=fc, ec=ec  )
        plt.fig.gca().add_artist(circle)

def plotGeom( atoms=None, atomSize=0.1 ):
    if atoms is not None:
        plotAtoms( atoms, atomSize=atomSize )

# --- plotting part here, plots all calculated signals --- #

NoV = len(didv) if didv_b else len(current)
NoH = len(didv[0]) if didv_b else len(current[0])

#print "DEBUG: Voltages", Voltages
#print "DEBUG: namez", namez
#print "DEBUG: NoV", NoV
#print "DEBUG: NoH", NoH

if PNG :
    print("We go to plotting ")
    for vv in range(NoV):
        for k in range(NoH):
            #print "DEBUG: long name:::", namez[vv],';height:%03d;tip:'  %k,tip_type,';',tip_orb
            name_plot=namez[vv]+';height:'+str(k)+';tip:'+tip_type+';'+tip_orb1+tip_orb2
            if didv_b :
                # ploting part here:
                plt.figure( figsize=(0.5 * lvec[1,0] , 0.5 * lvec[2,1] ) )
                plt.imshow(didv[vv,k,:,:], origin='image', extent=extent , cmap='gray')
                plotGeom(atoms=geom_plot)
                plt.xlabel(r' Tip_x $\AA$')
                plt.ylabel(r' Tip_y $\AA$')
                plt.title("dIdV:"+name_plot)
                plt.savefig( 'didv_'+namez[vv]+"_tip_"+tip_type+"-"+str(tip_orb1_amount)+tip_orb1+"-"+str(tip_orb2_amount)+tip_orb2+"_WF_"+str(WorkFunction+Voltages[vv]*WF_decay)+"_eta_"+str(eta)+'_%03d.png' %k , bbox_inches='tight' )
                plt.close()
            if STM_b :
                # ploting part here:
                plt.figure( figsize=(0.5 * lvec[1,0] , 0.5 * lvec[2,1] ) )
                plt.imshow(current[vv,k,:,:], origin='image', extent=extent , cmap='gray')
                plotGeom(atoms=geom_plot)
                plt.xlabel(r' Tip_x $\AA$')
                plt.ylabel(r' Tip_y $\AA$')
                plt.title("STM:"+name_plot)
                plt.savefig( 'STM_'+namez[vv]+"_tip_"+tip_type+"-"+str(tip_orb1_amount)+tip_orb1+"-"+str(tip_orb2_amount)+tip_orb2+"_WF_"+str(WorkFunction)+"_WF_decay_"+str(round(WF_decay,1))+"_eta_"+str(eta)+'_%03d.png' %k , bbox_inches='tight' )
                plt.close()
    print("Everything plotted")
if WSxM :
    print("writing WSxM files")
    for vv in range(NoV):
        for k in range(NoH):
            if didv_b :
                name_file =  'didv_'+namez[vv]+"_tip_"+tip_type+"-"+str(tip_orb1_amount)+tip_orb1+"-"+str(tip_orb2_amount)+tip_orb2+"_WF_"+str(WorkFunction+Voltages[vv]*WF_decay)+"_eta_"+str(eta)+'_%03d.xyz' %k
                tmp_curr=didv[vv,k,:,:].flatten()
                out_curr=np.zeros((len(tmp_curr),3))
                out_curr[:,0]=tip_r0[k,:,:,0].flatten()
                out_curr[:,1]=tip_r0[k,:,:,1].flatten()
                out_curr[:,2]=tmp_curr.copy()
                f=open(name_file,'w')
                print("WSxM file copyright Nanotec Electronica", file=f)
                print("WSxM ASCII XYZ file; obtained from dIdV code by Krejci et al.", file=f)
                print("X[A]  Y[A]  Z[A]", file=f)
                print("", file=f)
                np.savetxt(f, out_curr)
                f.close()
                #
            if STM_b :
                name_file =  'STM_'+namez[vv]+"_tip_"+tip_type+"-"+str(tip_orb1_amount)+tip_orb1+"-"+str(tip_orb2_amount)+tip_orb2+"_WF_"+str(WorkFunction)+"_WF_decay_"+str(round(WF_decay,1))+"_eta_"+str(eta)+'_%03d.xyz' %k
                tmp_curr=current[vv,k,:,:].flatten()
                out_curr=np.zeros((len(tmp_curr),3))
                out_curr[:,0]=tip_r0[k,:,:,0].flatten()
                out_curr[:,1]=tip_r0[k,:,:,1].flatten()
                out_curr[:,2]=tmp_curr.copy()
                f=open(name_file,'w')
                print("WSxM file copyright Nanotec Electronica", file=f)
                print("WSxM ASCII XYZ file; obtained from dIdV code by Krejci et al.", file=f)
                print("X[A]  Y[A]  Z[A]", file=f)
                print("", file=f)
                np.savetxt(f, out_curr)
                f.close()
                #
    print("WSxM files written")

if XSF :
    print("writing XSF files")
    xsf_head = Bu.At2XSF(geom_plot) if plot_atoms else GU.XSF_HEAD_DEFAULT
    for vv in range(NoV):
        if didv_b :
            name_file =  'didv_'+namez[vv]+"_tip_"+tip_type+"-"+str(tip_orb1_amount)+tip_orb1+"-"+str(tip_orb2_amount)+tip_orb2+"_WF_"+str(WorkFunction+Voltages[vv]*WF_decay)+"_eta_"+str(eta)+'.xsf'
            GU.saveXSF(name_file, didv[vv], lvec, head=xsf_head )
        if STM_b :
            name_file =  'STM_'+namez[vv]+"_tip_"+tip_type+"-"+str(tip_orb1_amount)+tip_orb1+"-"+str(tip_orb2_amount)+tip_orb2+"_WF_"+str(WorkFunction)+"_WF_decay_"+str(round(WF_decay,1))+"_eta_"+str(eta)+'.xsf'
            GU.saveXSF(name_file, current[vv], lvec, head=xsf_head )
    print("XSF files written")

if NPY :
    print("writing npy binary files")
    for vv in range(NoV):
        if didv_b :
            name_file =  'didv_'+namez[vv]+"_tip_"+tip_type+"-"+str(tip_orb1_amount)+tip_orb1+"-"+str(tip_orb2_amount)+tip_orb2+"_WF_"+str(WorkFunction+Voltages[vv]*WF_decay)+"_eta_"+str(eta)
            GU.saveNpy(name_file, didv[vv], lvec)#, head=XSF_HEAD_DEFAULT )
        if STM_b :
            name_file =  'STM_'+namez[vv]+"_tip_"+tip_type+"-"+str(tip_orb1_amount)+tip_orb1+"-"+str(tip_orb2_amount)+tip_orb2+"_WF_"+str(WorkFunction)+"_WF_decay_"+str(round(WF_decay,1))+"_eta_"+str(eta)
            GU.saveNpy(name_file, current[vv], lvec)#, head=XSF_HEAD_DEFAULT )
    print("npy files written")

# --- the end --- #

print() 
print()
print("Done")
print()
