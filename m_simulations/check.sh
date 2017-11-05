#!/bin/bash
for i in `seq 1 500`
do
if [ -f "./results/GP/$@/GP_$@_data_m_$i.RData" ]
then
echo "GP_$@_data_m_$i.RData is found."
else
qsub -v casenumber=`printf %03d $i`,distr=$@ -N $@$i case.pbs;
fi
done
