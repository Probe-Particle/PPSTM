#!/bin/bash

# Tested only on Linux
# !!! git & internet needed !!!
# This script will download for you ProbeParticleModel,
# it is necessary for smulations with tilting tip

directory=${PWD##*/}
cd ..
mkdir PPAFM
cd PPAFM
git init
git pull https://github.com/konradkusiak97/Aalto-SIN-PPAFM-python3.git
cd ../$directory
ln -s ../PPAFM ./PPAFM
ln -s ../PPAFM/pyProbeParticle pyProbeParticle