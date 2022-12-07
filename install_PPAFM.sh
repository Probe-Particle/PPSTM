#!/bin/bash

# Tested only on Linux
# !!! git & internet needed !!!
# This script will download for you ProbeParticleModel,
# it is necessary for smulations with tilting tip
# This script downloads the old master branch of the PPAFM code

git clone -b master_backup git@github.com:Probe-Particle/ProbeParticleModel.git ../PPAFM 
ln -s ../PPAFM ./PPAFM
ln -s ../PPAFM/pyProbeParticle ./pyProbeParticle