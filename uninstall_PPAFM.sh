#!/bin/bash

# Tested only on Linux
# !!! git & internet needed !!!
# This script will download for you ProbeParticleModel,
# it is necessary for smulations with tilting tip

directory=${PWD##*/}
rm -r PPAFM
rm -r pyProbeParticle
cd ..
rm -r PPAFM