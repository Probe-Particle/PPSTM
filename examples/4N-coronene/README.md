# HOW TO RUN PPSTM simulations:
Relaxed scans require installation of PPAFM (version =< 0.2.0a3) with 
`pip install ppafm==0.2.0a3` # See the [ppafm installation](https://github.com/Probe-Particle/ppafm/wiki/Install-ppafm) if you run into problems 

```
ppafm-generate-ljff -i crazy_mol.xyz -f npy # see params.ini &/or AFM documentation for initial settings.
```
And now LJ Force-field is made. (This simple examplary run does not use charget tip-apex).

```
ppafm-relaxed-scan --pos -f npy # AFM scan with CO tip, positions of PP saved in Q*.**K*.**/PPpos_?.npy
ppafm-plot-results --pos --df --save_df -f npy # save figures of PP positions; figures of df and save df into Q*.**K*.**/Amp*.**/df.pny file 
```

Now all AFM files are made; in the minimalistic mode only first two commands are needed. But it is always good to check, whether one is
in proper height, by means of looking at xy and df, before running dI/dV or even STM, which in case of slab calculations can last long time
Before running the script, please look at the dIdV script. For more questions look at the [wiki documentation](https://github.com/Probe-Particle/PPSTM/wiki)
or [contact us](https://github.com/Probe-Particle/PPSTM/wiki#contact). Now running dI/dV scan over HOMO and LUMO on the electronic structure
precalculated by Fireball.:

```
python PPSTM/dIdV_test_4N-coronene.py # This is an old script which will probably be shortly changed with ppstm_run ****.toml
```

Simulation over the molekule at energies HOMO and LUMO are plotted. 1st figure of four is AFM, 2nd STM without relaxation with s-tip:
that means it reveals the original orbital; 3rd relaxed STM scan with s-orbital on PP; 4th relaxed STM with px & py orbital on PP.
Please note that df is done with oscilating tip, whether STM are not. The height between AFM and STM is therefore aproximative in the figures.
To calculate the averaged current is in principle possible, but to show RAW data from single height is more straightforward according to us.

## Precalculations by GPAW (needs ASE & GPAW installed):
```
python run_GPAW_LCAO.py  # .gpw file created - basicly what is needed for the STM calculations
python print_wf_sp.py    # creates phik_*.dat with written LCAO coefficients for the simulations (optional, nod needed, since PPSTM can read GPAW
			 # inputs (*.gpw) within ppstm_run.py automatically
python plot_frontiers.py # print cube files with HOMO & LUMO wf for visualization
```

## Precalculations by CP2K (needs a CP2K installed and additional files with basis and potentials):
```
cp2k crazy_mol_cp2k.inp  # the input file was created around 2016, so we cannot ensure that it works with the latest version of cp2k
```
This will create `crazy_mol-cartesian-mos-1_0.MOLog` which have the electronic structure written in. After that you can use:

```
python PPSTM/dIdV_test_4N-coronene_cp2k.py # This is an old script which will probably be shortly changed with ppstm_run ****.toml
```
