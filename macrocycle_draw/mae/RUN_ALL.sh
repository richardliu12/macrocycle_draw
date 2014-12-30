#!/bin/bash

# submits conf. search jobs as job array
# saves in ./output/$filename$/

COUNTER=0

for i in *.com; do
    let "COUNTER++" 
    cp $i csearch"$COUNTER".com
done

cp ../slurm_array.sh ./slurm_array.sh
if ($COUNTER > 4)
then
    sbatch --array=1-"$COUNTER"\%4 slurm_array.sh
else
    sbatch slurm_array.sh
fi
rm -f slurm_array.sh
