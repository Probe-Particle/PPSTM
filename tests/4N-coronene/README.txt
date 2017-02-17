#HOW TO RUN PPSTM simulations:
#See doc, or:


python PPSTM/PPAFM/generateLJFF.py -i crazy_mol.xyz --npy # see params.ini &/or AFM documentation for initial settings
# now LJ Force-field is made

python PPSTM/PPAFM/relaxed_scan.py --pos --npy # AFM scan with CO tip, positions of PP saved in Q*.**K*.**/PPpos_?.npy
python PPSTM/PPAFM/plot_results.py --pos --df --save_df --npy # save figures of PP positions; figures of df and save df into Q*.**K*.**/Amp*.**/df.pny file 

# Now all AFM files are made; in the minimalistic mode only first two commands are needed. But it is always good to check, whether one is
# in proper height, by means of looking at xy and df, before running dI/dV or even STM, which in case of slab calculations can last long time
# Before running the script, please look at the dIdV script. For more questions look at the documentation or contact me at krejcio@fzu.cz
# running dI/dV scan over HOMO and LUMO:

python PPSTM/dIdV_test_4N-coronene.py

# Simulation over the molekule at energies HOMO and LUMO are plotted. 1st figure of four is AFM, 2nd STM without relaxation with s-tip:
# that means it reveals the original orbital; 3rd relaxed STM scan with s-orbital on PP; 4th relaxed STM with px & py orbital on PP.
# Please note that df is done with oscilating tip, whether STM are not. The height between AFM and STM is therefore aproximative in the figures.
# To calculate the averaged current is in principle possible, but to show RAW data from single height is more straightforward according to us.

## Precalculations by GPAW (needs ASE & GPAW installed):
python run_GPAW_LCAO.py  # .gpw file created - basicly what is needed for the STM calculations
python print_wf_sp.py    # creates phik_*.dat with written LCAO coefficients for the simulations (optional, nod needed, since PPSTM can read GPAW
			 # imputs imidiatelly
python plot_frontiers.py # print cube files with HOMO & LUMO wf


