#!/bin/bash
set -euo pipefail

# resetting plot_file
> ATLAS_2020_I1803608.plot
cat ATLAS_2020_I1803608_original.plot > ATLAS_2020_I1803608.plot

makeplots0="BEGIN PLOT /ATLAS_2020_I1803608:CUTS=$1:TYPE=EW_ONLY/"


# adding plots
echo "" >> ATLAS_2020_I1803608.plot

for i in $(seq 0 7); do
cat >> ATLAS_2020_I1803608.plot << EOF
${makeplots0}A_$i
XLABEL=\$p_T\$ [GeV]
YLABEL=\$A_$i\$
LogY=0
END PLOT

EOF
done
