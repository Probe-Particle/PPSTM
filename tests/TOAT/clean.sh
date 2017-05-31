#!/bin/bash

echo "Removing not neccessary folders and *.npy files"
rm -r ./Q0.00K0.50
rm *.npy
echo "Removing calculated Figures"
rm didv_*.png
echo "Done, bye bye"


