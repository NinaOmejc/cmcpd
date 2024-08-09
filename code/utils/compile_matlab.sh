#!/bin/sh
#
# Compile the matlab code so we can run it without a matlab license. Required:
#     Matlab 2023a including compiler, with license

export PATH=/usr/local/MATLAB/R2023a/bin:$PATH

mcc -m -v -I /home/nina/Experiments/corticomuscular_analysis/src/AMICA1.7 src/main.m -d bin
chmod go+rx bin/main
chmod go+rx bin/run_main.sh


