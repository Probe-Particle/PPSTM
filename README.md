# PPSTM (Probe Particle STM)
Code simulating various various STM techniques, especially for [tilting tips](https://pubs.acs.org/doi/10.1021/ja204624g) (depending on [ProkopHapala/ProbeParticleModel](https://github.com/ProkopHapala/ProbeParticleModel) )
This is implementation of efficient and simple model for simulation of High-resolution scanning tunneling microscopy (STM).
Normall STM simulations using [Chen's approximattion](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.42.8841) is part of the code as well

- BEWARE - after repairing misteake in d-orbitals the results for FePc are different. We are investigating it.

* currently developed Python3/C++ version of the PPSTM code:
  * works in python3(.7 and higher), matplotlib, cpp>=4.4.8 (look at wiki instruction for problems with MAC [compilation](https://github.com/Probe-Particle/PPSTM/wiki#compilation-and-overview)), PyQt5 (for [GUI](https://github.com/Probe-Particle/PPSTM/wiki#GUI-for-PPSTM-code) ),for some parts [ASE](https://wiki.fysik.dtu.dk/ase/) and [GPAW](https://wiki.fysik.dtu.dk/gpaw/) are imporatant;
  * Part of the code regarding simulations with tilting tips is depending on the ProbeParticleModel developed by Prokop Hapala and co. (https://github.com/ProkopHapala/ProbeParticleModel), you can easilly install it throug ```install_PPAFM.sh``` script.

For easy introduction to this code and its functionalities, try [Graphic User Interface (GUI)](https://github.com/Probe-Particle/PPSTM/wiki#GUI-for-PPSTM-code) or [PPSTM_simple.py](https://github.com/Probe-Particle/PPSTM/wiki#ppstm_simplepy) script.

**Documentation is here at [Wiki](https://github.com/Probe-Particle/PPSTM/wiki).**

It can also simulate IETS images of molecules, if the imaging mechanism is driven by the amplitude of the IETS peak.

### References (STM)
* [Ondrej Krejčí, Prokop Hapala, Martin Ondráček, and Pavel Jelínek, Principles and simulations of high-resolution STM imaging with a flexible tip apex, Phys. Rev. B 95, 045407 – Published 6 January 2017 ](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.95.045407) 

### References (IETS & STM-d orbitals)
* [Bruno de la Torre, Martin Švec, Giuseppe Foti, Ondřej Krejčí, Prokop Hapala, Aran Garcia-Lekue, Thomas Frederiksen, Radek Zbořil, Andres Arnau, Héctor Vázquez, and Pavel Jelínek, Submolecular Resolution by Variation of the Inelastic Electron Tunneling Spectroscopy Amplitude and its Relation to the AFM/STM Signal, Phys. Rev. Lett. 119, 166001 – Published 16 October 2017](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.119.166001)
