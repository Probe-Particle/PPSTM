#!/bin/bash

echo "Removing calculated Figures"
rm PDOS_*.png
rm didV-*.png
rm STM-*.png
rm STM_*.png
echo "Removing outputs of PPSTM_simple, if any"
rm STM_*.png
rm didv_*.png
rm didv_*tip_*WF_*eta_*.xsf
rm STM_*tip_*WF_*eta_*.xsf
rm didv_*tip_*WF_*eta_*.xyz
rm STM_*tip_*WF_*eta_*.xyz
rm PDOS_*.png
echo "Done, bye bye"


