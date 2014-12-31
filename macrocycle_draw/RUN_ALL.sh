#!/bin/bash

# submits conf. search jobs as job array
# saves in ./output/$filename$/

maxjobs=4

for i in *.com; do
    i=`echo $i | cut -f1 -d '.'`
    ./RUN_ONE.sh ${i} &
    echo Running $i
    while [[ $(jobs|wc -l) -ge "$maxjobs" ]]
        do sleep 10
    done
done
