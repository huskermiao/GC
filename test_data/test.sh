#!/bin/sh

set -x

for i in F2_120i_6chr_0.03dens.map F2_X0.45miss0.35.map RIL_141i_8chr.map RIL_X0.4miss0.35.map
do
    perl ../binmarkers.pl -e X -o A,B -t 5 $i
done
