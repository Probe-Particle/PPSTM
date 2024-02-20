import os
import pdb
import sys
import time

import numpy as np
import h5py

import ppafm

class STMgenerator():
    """
    Class for creating STM input/output pairs for machine learning.
    An Iterator.

    This class takes as parameters the location of the input h5 file,
    specifications for the data set,(e.g. batch size and the output type) 
    and an instance of the STMulator.STMulator class.

    Arguments:
        stmulator: STMulator.STMulator.
        auxmaps: List of AuxMap objects.
        h5_name: string. Path to the h5 file containing the molecule data
        rotations_pkl: A pickled npy file containing a dictionary with rotations for each molecule. {mode: {cid: [rotations]}}
        mode: string. train/val/test. Which data set to create. (default='train')
        batch_size: int. Batch size for the data set. (default=32)
        dist_above: float. The tip-sample distance at the lowest tip position. (default=5.3)
        timings: boolean. Prints timing information if set True. (default=False)
    """

    def __init__(self, 
        stmulator, auxmaps, h5_name, rotations_pkl,
        mode = 'train',
        batch_size = 32,
        dist_above = 5.0,
        timings=False
        ):

        self.stmulator = stmulator
        self.auxmaps = auxmaps
        self.h5 = h5_name
        self.rotations_pkl = rotations_pkl
        self.mode = mode
        self.batch_size = batch_size
        self.dist_above = dist_above
        self.dist_above_active = dist_above

        self.rotations = []
        self.read_rotations_and_shuffle()
        self.counter = 0
        self.timings = timings

    def __next__(self):
        if self.counter < len(self.rotations):
            molecules = []
            X1s = []
            X2s = []

            Ys = [[] for _ in range(len(self.auxmaps))]
            batch_size = min(self.batch_size, len(self.rotations) - self.counter)

            if self.timings: batch_start = time.time()

            for i in range(batch_size):
                if self.timings: sample_start = time.time()
                
                # Select cid + rotation with index self.counter
                cid, rot = self.rotations[self.counter]

                # Load data from h5 file
                if self.timings: h5_start = time.time()
                eigs, coefs, xyz, qs, Zs, lvec = self.get_data_from_h5(cid)
                if self.timings: print(f"Data loading time: {time.time()-h5_start:.2f} s")

                # Rotate xyz
                xyz = self.rotate_molecule(xyz, rot)

                self.xyz = xyz
                self.qs = qs
                self.Zs = Zs

                # Handle position + distance
                self.handle_positions()
                self.handle_distance()

                # Run afm + stm, tip relaxing done internally
                x_stm, x_afm = self.stmulator(self.xyz, self.Zs, self.qs, eigs, coefs, REAs=self.REAs)
                X1s.append([x_stm])
                X2s.append([x_afm])

                self.xyz[:, [1,0]] = self.xyz[:, [0,1]]

                mol = np.c_[self.xyz, self.Zs, self.qs]
                molecules.append(mol)
                
                # Calculate each aux map
                for ia, auxmap in enumerate(self.auxmaps):
                    xyzqs = np.concatenate([self.xyz, self.qs[:,None]], axis=1)
                    Ys[ia].append(auxmap(xyzqs, self.Zs))

                if self.timings: print(f'Sample {i} runtime [s]: {time.time() - sample_start}')
                self.counter += 1
            
            for i in range(len(self.auxmaps)):
                Ys[i] = np.stack(Ys[i], axis=0)

            if self.timings: print(f'Batch runtime [s]: {time.time() - batch_start}')

        else:
            raise StopIteration

        return (X1s, X2s), Ys,  molecules

    
    def __iter__(self):
        self.counter = 0
        return self

    def __len__(self):
        """Returns the number of batches to be generated."""
        return len(self.rotations)//self.batch_size

    def read_rotations_and_shuffle(self):
        """Reads the file containing rotations into a list of tuples and shuffles it."""
        d = np.load(self.rotations_pkl, allow_pickle=True)
        cid_rot_list = [(cid, rot) for cid in d[self.mode].keys() for rot in d[self.mode][cid]]
        r = np.random.RandomState(0) # use a predefined random state for reproduciblity
        r.shuffle(cid_rot_list)
        self.rotations = cid_rot_list


    def get_data_from_h5(self, cid):
        """
        Reads properties of a given molecule in the h5 file.

        Arguments:
            cid: string or int. The identification number for the molecule.
        Returns: 6-tuple. Eigenenergies, KS eigenvectors, atomic positions, atomic charges, atomic numbers and the unit cell boundaries.
        """
        if isinstance(cid, int): cid = str(cid)

        with h5py.File(self.h5, 'r') as f:
            g = f[cid]
            eigs = np.array(g['eigs'])
            coefs = np.array(g['coefs'])
            xyz = g['xyz']
            lattice = np.array(xyz.attrs['lattice'])
            xyz = np.array(xyz)
            Zs = np.array(g['Z'])
            qs = np.array(g['qs'])
        return eigs, coefs, xyz, qs, Zs, lattice

    def handle_positions(self):
        '''
        Set current molecule to the center of the scan window.
        '''
        sw = self.stmulator.scan_window
        scan_center = np.array([sw[1][0] + sw[0][0], sw[1][1] + sw[0][1]]) / 2
        self.xyz[:,:2] += scan_center - self.xyz[:,:2].mean(axis=0)

    def handle_distance(self):
        '''
        Set correct distance from scan region for the current molecule.
        '''
        RvdwPP = self.stmulator.afmulator.typeParams[self.stmulator.afmulator.iZPP-1][0]
        self.REAs = ppafm.getAtomsREA(self.stmulator.iZPP, self.Zs, self.stmulator.afmulator.typeParams, alphaFac=-1.0)
        Rvdw = self.REAs[:,0] - RvdwPP
        zs = self.xyz[:,2]
        imax = np.argmax(zs + Rvdw)
        total_distance = self.dist_above_active + Rvdw[imax] + RvdwPP - (zs.max() - zs[imax])
        self.xyz[:,2] += (self.stmulator.scan_window[1][2] - total_distance) - zs.max()

    
    def rotate_molecule(self, xyz, rotation_matrix):
        """Rotates the atomic coordinates according to the rotation matrix."""
        rotated = np.dot(xyz, rotation_matrix.T)
        return rotated








