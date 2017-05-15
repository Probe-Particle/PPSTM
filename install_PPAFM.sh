#!\bin\bash

# Tested only on Linux
# !!! git & internet needed !!!
# This script will download for you ProbeParticleModel,
# it is necessary for smulations with tilting tip

cd ..
git pull https://github.com/ProkopHapala/ProbeParticleModel/
mv ProbeParticleModel  PPAFM
ln -s PPAFM PPSTM/PPAFM
cd PPSTM
ln -s PPAFM/pyProbeParticle pyProbeParticle