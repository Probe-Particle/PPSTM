#!/bin/bash

echo "Removing not neccessary folders and *.npy files"
rm -r ./Q0.00K0.50
rm *.npy
echo "Removing calculated Figures"
rm didv_*.png
echo "Removing outputs of PPSTM_simple, if any"
rm STM_*.png
rm didv_*tip_*WF_*eta_*.xsf
rm STM_*tip_*WF_*eta_*.xsf
rm didv_*tip_*WF_*eta_*.xyz
rm STM_*tip_*WF_*eta_*.xyz
echo "Done, bye bye"


