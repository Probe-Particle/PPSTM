[![Ask DeepWiki](https://deepwiki.com/badge.svg)](https://deepwiki.com/Probe-Particle/PPSTM)

# PPSTM (Probe Particle STM)
Code simulating various STM techniques, especially for [tilting tips](https://pubs.acs.org/doi/10.1021/ja204624g) (depending on [ProkopHapala/ProbeParticleModel](https://github.com/ProkopHapala/ProbeParticleModel) )
This is implementation of efficient and simple model for simulation of High-resolution scanning tunneling microscopy (STM).
Normall STM simulations using [Chen's approximattion](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.42.8841) is part of the code as well

- BEWARE - after repairing misteake in d-orbitals the results for FePc are different. We are investigating it.

* currently developed Python3/C++ version of the PPSTM code:
  * works in python3.12, matplotlib, cpp>=4.4.8 (look at wiki instruction for problems with [MAC compilation](https://github.com/Probe-Particle/PPSTM/wiki/Installation#mac)), PyQt5 (for [GUI](https://github.com/Probe-Particle/PPSTM/wiki#GUI-for-PPSTM-code) ),for some parts [ASE](https://wiki.fysik.dtu.dk/ase/) and [GPAW](https://wiki.fysik.dtu.dk/gpaw/) are imporatant;
  * Part of the code regarding simulations with tilting tips is depending on the PPAFM developed by Prokop Hapala and co. (https://github.com/Probe-Particle/ppafm), you can easilly install it with ```pip install ppafm>=0.2.0a3```.

For easy introduction to this code and its functionalities, try [Graphic User Interface (GUI)](https://github.com/Probe-Particle/PPSTM/wiki#GUI-for-PPSTM-code) or [PPSTM_simple.py](https://github.com/Probe-Particle/PPSTM/wiki#ppstm_simplepy) script.

Constant current simulations are now available through `const_cur_tutorial.ipynb` jupyter notebook from npz or xsf. There is also former [Mathematica notebook](https://github.com/Probe-Particle/MathematicaForPPSTM/blob/master/PPSTM_contant_current_XSF_view.nb) which works only with xsf.

### Documentation

**Documentation is here at [Wiki](https://github.com/Probe-Particle/PPSTM/wiki).**

It can also simulate IETS images of molecules, if the imaging mechanism is driven by the amplitude of the IETS peak.

#### Auto-Generated Documentation

* [![Ask DeepWiki](https://deepwiki.com/badge.svg)](https://deepwiki.com/Probe-Particle/PPSTM)

### Setup
A simple setup with CONDA [Conda](https://conda.io) on Linux is:.
```bash
conda env create -f environment.yml
conda activate ppstm
python -m ipykernel install --user --name=ppstm
```

More details can be found on a [dedicated Wikipage](https://github.com/Probe-Particle/PPSTM/wiki/Installation).

### Tests
After activating the `ppstm-dev` environment, you can run the test suite with:
```bash
pytest tests -v
```
This runs all tests in verbose mode.

#### Test Coverage
Measure code coverage to identify untested code paths.

##### View coverage in the terminal
To show coverage statistics by file and list lines not covered by tests:
```bash
pytest --cov=pyPPSTM --cov-report term-missing tests
```

##### Generate an interactive HTML report
To create a detailed HTML report with coverage by file, function, and class:
```bash
pytest --cov=pyPPSTM --cov-report html:coverage tests
```
Open `coverage/index.html` in your browser to explore the results interactively.

### References (should be always cited)
* [Ondrej Krejčí, Prokop Hapala, Martin Ondráček, and Pavel Jelínek, Principles and simulations of high-resolution STM imaging with a flexible tip apex, Phys. Rev. B 95, 045407 – Published 6 January 2017 ](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.95.045407) 

### References (should be cited for the relaxed scans)
* [Niko Oinonen, Aliaksandr V. Yakutovich, Aurelio Gallardo, Martin Ondracek, Prokop Hapala, Ondrej Krejci, Advancing Scanning Probe Microscopy Simulations: A Decade of Development in Probe-Particle Models, Comput. Phys. Commun. 305, 109341 - Available online 10 August 2024](https://doi.org/10.1016/j.cpc.2024.109341)

### References (should be cited if IETS or STM-d orbitals were included) 
* [Bruno de la Torre, Martin Švec, Giuseppe Foti, Ondřej Krejčí, Prokop Hapala, Aran Garcia-Lekue, Thomas Frederiksen, Radek Zbořil, Andres Arnau, Héctor Vázquez, and Pavel Jelínek, Submolecular Resolution by Variation of the Inelastic Electron Tunneling Spectroscopy Amplitude and its Relation to the AFM/STM Signal, Phys. Rev. Lett. 119, 166001 – Published 16 October 2017](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.119.166001)

