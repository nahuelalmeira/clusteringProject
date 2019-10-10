#!/bin/bash

#$ -cwd
#$ -l h_rt=120:00:00


N=4000
alpha="020"
transit=100
decorr=10

/bin/hostname
/bin/date

./HM_cycle.sh ${N} ${alpha} ${transit} ${decorr} 

echo " "
echo "Job finished at: "
/bin/date
