import numpy as np
import matplotlib.pyplot as plt

from pyPPSTM import elements
from pyPPSTM import basUtils as Bu


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

def get_voltages_and_names(config, eigs):
    scan_type = config['scan']['scan_type']

    V = config['scan']['V']
    V_max = config['scan']['V_max']
    dV = config['scan']['dV']

    if scan_type in ['states', 'STATES']:
        states = np.sort(eigs)
        mask = (states >= V) & (states <= V_max)
        voltages = states[mask]
        names = [f"{V:.5f}" for V in voltages]
    else:
        voltages = np.arange(V, V_max+0.001 ,dV)
        names = [f"{V:.1f}" for V in voltages]

    return voltages, names

def get_number_of_voltages_and_heights(config, current, didv):
    scan_type = config['scan']['scan_type']

    if scan_type in ['STM', 'STM-single']:
        nV = current.shape[0]
        nH = current.shape[1]
    else:
        nV = didv.shape[0]
        nH = didv.shape[1]

    return nV, nH

def plot_png(config, current, didv, voltages, names, lvec, extent, geom_plot):
    tip_type = config['scan']['tip_type']
    tip_orb = config['scan']['tip_orb']
    eta = config['scan']['eta']
    work_function = config['advanced']['work_function']
    wf_decay = config['advanced']['work_function_decay']

    nV, nH = get_number_of_voltages_and_heights(config, current, didv)
    for vv in range(nV):
        for k in range(nH):
            #print "DEBUG: long name:::", namez[vv],';height:%03d;tip:'  %k,tip_type,';',tip_orb
            name_plot = f'{names[vv]};height:{k};tip:{tip_type};{tip_orb}'
            if didv is not None:
                # ploting part here:
                plt.figure(figsize=(0.5*lvec[1,0], 0.5*lvec[2,1]))
                plt.imshow(didv[vv,k,:,:], origin='lower', extent=extent, cmap='gray')
                plotGeom(atoms=geom_plot)
                plt.xlabel(r' Tip_x $\AA$')
                plt.ylabel(r' Tip_y $\AA$')
                plt.title("dIdV:"+name_plot)
                save_name = f'didv_{names[vv]}_tip_{tip_type}-{tip_orb}_WF_{work_function-voltages[vv]*wf_decay}_eta_{eta:.1f}_{k:03d}.png'
                plt.savefig(save_name, bbox_inches='tight')
                plt.close()
            if current is not None:
                # ploting part here:
                plt.figure(figsize=(0.5*lvec[1,0], 0.5*lvec[2,1]))
                plt.imshow(current[vv,k,:,:], origin='lower', extent=extent, cmap='gray')
                plotGeom(atoms=geom_plot)
                plt.xlabel(r' Tip_x $\AA$')
                plt.ylabel(r' Tip_y $\AA$')
                plt.title("STM:"+name_plot)
                save_name = f'STM_{names[vv]}_tip_{tip_type}-{tip_orb}_WF_{work_function:.1f}_WF_decay_{wf_decay:.1f}_eta_{eta:.1f}_{k:03d}.png'
                plt.savefig(save_name, bbox_inches='tight')
                plt.close()

def plot_wsxm(config, current, didv, voltages, names, tip_r0):
    tip_type = config['scan']['tip_type']
    tip_orb = config['scan']['tip_orb']
    eta = config['scan']['eta']
    work_function = config['advanced']['work_function']
    wf_decay = config['advanced']['work_function_decay']
    nV, nH = get_number_of_voltages_and_heights(config, current, didv)

    for vv in range(nV):
        for k in range(nH):
            if didv is not None:
                name_file = f'didv_{names[vv]}_tip_{tip_type}-{tip_orb}_WF_{work_function-voltages[vv]*wf_decay}_eta_{eta:.1f}_{k:03d}.wsxm'
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
        
            if current is not None:
                name_file = f'STM_{names[vv]}_tip_{tip_type}-{tip_orb}_WF_{work_function:.1f}_WF_decay_{wf_decay:.1f}_eta_{eta:.1f}_{k:03d}.wsxm'
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

def save_xsf(config, current, didv, voltages, names, geom_plot, lvec):
    try:
        import ppafm.io as io
    except ImportError:
        raise ImportError("ppafm needed for XSF output")
    
    tip_type = config['scan']['tip_type']
    tip_orb = config['scan']['tip_orb']
    eta = config['scan']['eta']
    work_function = config['advanced']['work_function']
    wf_decay = config['advanced']['work_function_decay']

    xsf_head = Bu.At2XSF(geom_plot) if geom_plot is not None else io.XSF_HEAD_DEFAULT
    nV, nH = get_number_of_voltages_and_heights(config, current, didv)
    
    for vv in range(nV):
        if didv is not None:
            name_file = f'didv_{names[vv]}_tip_{tip_type}-{tip_orb}_WF_{work_function-voltages[vv]*wf_decay}_eta_{eta:.1f}.xsf'
            io.saveXSF(name_file, didv[vv], lvec, head=xsf_head )
        if current is not None:
            name_file = f'STM_{names[vv]}_tip_{tip_type}-{tip_orb}_WF_{work_function:.1f}_WF_decay_{wf_decay:.1f}_eta_{eta:.1f}.xsf'
            io.saveXSF(name_file, current[vv], lvec, head=xsf_head )

def save_npy(config, current, didv, voltages, names, lvec, atomic_info_or_head):
    try:
        import ppafm.io as io
    except ImportError:
        raise ImportError("ppafm needed for NPY output")
    tip_type = config['scan']['tip_type']
    tip_orb = config['scan']['tip_orb']
    eta = config['scan']['eta']
    work_function = config['advanced']['work_function']
    wf_decay = config['advanced']['work_function_decay']

    nV, nH = get_number_of_voltages_and_heights(config, current, didv)

    for vv in range(nV):
        if didv is not None:
            name_file = f'didv_{names[vv]}_tip_{tip_type}-{tip_orb}_WF_{work_function-voltages[vv]*wf_decay}_eta_{eta:.1f}'
            io.saveNpy(name_file, didv[vv], lvec, atomic_info=atomic_info_or_head)
        if current is not None:
            name_file = f'STM_{names[vv]}_tip_{tip_type}-{tip_orb}_WF_{work_function:.1f}_WF_decay_{wf_decay:.1f}_eta_{eta:.1f}'
            io.saveNpy(name_file, current[vv], lvec, atomic_info=atomic_info_or_head)
