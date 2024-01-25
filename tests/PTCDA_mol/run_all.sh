#!/bin/bash

########################################################################
#                                                                      #
# This is a script that will create the PPSTM simulations of HOMO and  #
# LUMO of PTCDA molecule (in vacuum) as published in the paper.        #
# First it will run the PPAFM to create the positions of probe particle#
# that are then used for calculation of the dIdV masp                  #
#                                                                      #
########################################################################

wget https://zenodo.org/records/10563098/files/Fig_4_STM_data.tar.gz
tar -xzvf Fig_4_STM_data.tar.gz
ppafm-generate-ljff -i cube_001_hartree_potential.cube
ppafm-generate-elff -i cube_001_hartree_potential.cube
ppafm-relaxed-scan --pos
python3 xsf2xyz.py
sed -i '2d' input_plot.xyz # the xyz reader in PPSTM does not like non-empty 2nd line of the xyz file #
python3 PPSTM_homo.py
python3 PPSTM_lumo.py
