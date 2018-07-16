#!/bin/bash

p_lowPSD=0.19
d_lowPSD=0.19
p_lowE=0.57
d_lowE=0.66
zLow=-1000
zHigh=1000
timeBin=23.5
dtFit=1

root -l -b <<EOF 

.L Calculate/RnPoVsCell.C+
RnPoVsCell($p_lowPSD,$d_lowPSD,$p_lowE,$d_lowE,$zLow,$zHigh,$dtFit)

.q
EOF

root -l -b <<EOF 

.L Calculate/RnPoVsTime.C+
RnPoVsTime($p_lowPSD,$d_lowPSD,$p_lowE,$d_lowE,$zLow,$zHigh,$timeBin,$dtFit)

.q
EOF

#root -l -b <<EOF 

#.L Calculate/RnPoColVsTime.C+
#RnPoColVsTime($p_lowPSD,$d_lowPSD,$p_lowE,$d_lowE,$timeBin,$dtFit)

#.q
#EOF

#root -l -b <<EOF 

#.L Calculate/RnPoRowVsTime.C+
#RnPoRowVsTime($p_lowPSD,$d_lowPSD,$p_lowE,$d_lowE,$timeBin,$dtFit)

#.q
#EOF



