#!/bin/bash

# Tested only on Linux
# !!! git & internet needed !!!
# This script will download for you ProbeParticleModel,
# it is necessary for smulations with tilting tip

directory=${PWD##*/}
cd ..
mkdir PPAFM
cd PPAFM
#git init
git clone https://github.com/ProkopHapala/ProbeParticleModel
#git clone https://github.com/ProkopHapala/ProbeParticleModel -b OpenCL_py3 ./
cd ../$directory
ln -s ../PPAFM ./PPAFM
ln -s ../PPAFM/pyProbeParticle pyProbeParticle
