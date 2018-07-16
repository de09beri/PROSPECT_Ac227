#!/bin/bash

i=0;

while [ $i -lt 154  ]
do

root -l -b <<EOF
.L PlotDistributionsVsCell.C+
PlotDistributionsVsCell($i)
.q

EOF

i=$[$i+1]
done
