#!/bin/bash

echo "Removing calculated Figures"
rm didV-*.png
echo "Removing outputs of PPSTM_simple, if any"
rm STM_*.png
rm didv_*.png
rm didv_*tip_*WF_*eta_*.xsf
rm STM_*tip_*WF_*eta_*.xsf
rm didv_*tip_*WF_*eta_*.xyz
rm STM_*tip_*WF_*eta_*.xyz
echo "Done, bye bye"


