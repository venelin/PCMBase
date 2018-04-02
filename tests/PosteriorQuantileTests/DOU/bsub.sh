#!/bin/bash
bsub -W 23:59 -n 24 -R "rusage[mem=5000]" -J "LR_asim_rep" "R --vanilla --slave -f BayesValidateTwoSpeedOU.R"
