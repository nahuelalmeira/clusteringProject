#!/bin/bash

N=$1
alpha=$2
MODE="maxC"
randMCS=100
inputDir="../../../networks/model/HM/HM_m4/HM_m4_N${N}/HM_m4_N${N}_alpha${alpha}"
NETWORK="HM_m4_N${N}_alpha${alpha}_00000_gcc"
nSamples=100
TRANSIT=$3
DECORR=$4
SEED=0

../annealing ${NETWORK} --inputDir ${inputDir} --mode ${MODE} --randMCS ${randMCS} \
                        --transitTime ${TRANSIT} --decorrTime ${DECORR} \
                        --samples ${nSamples}  --cycle --seed ${SEED}
