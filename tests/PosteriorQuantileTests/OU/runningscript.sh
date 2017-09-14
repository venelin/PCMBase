#!/bin/bash
for i in {49..96}
do
bsub -W 23:59 -n 1 -R "rusage[mem=5000]" -J "LR_asim_rep" "R --vanilla --slave -f runReplication.R --args $i"
done
exit 0
