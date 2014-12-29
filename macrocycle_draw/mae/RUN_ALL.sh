#!/bin/bash

# submits conf. search jobs as job array
# saves in ./output/$filename$/

COUNTER=1

for i in *.com; do
    cp $i "csearch$COUNTER.com"
done

