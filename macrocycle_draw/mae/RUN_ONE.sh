#!/bin/bash

# $1 should not have a file extension
# submits conf. search $1.com
# saves in ./output/$1/

# run csearch
date
/n/sw/schrodinger/macromodel -LOCAL -WAIT $1.com &
wait

# set up minimization on results
cd ../output/$1
mv $1-csearch.maegz $1-csearch.gz
gunzip $1-csearch.gz
mv $1-csearch $1-csearch.mae

echo "Minimizing found conformations..."
inputName=$1"-csearch.mae"
outputName=$1"-csearch-min.maegz"
sed 1s/.*/"$inputName"/ ../../multimini.com | sed 2s/.*/"$outputName"/ > $1-csearch-min.com

# run minimizations
/n/sw/schrodinger/macromodel -LOCAL -WAIT $1-csearch-min.com
mv $1-csearch-min.maegz $1-csearch-min.gz
gunzip $1-csearch-min.gz
mv $1-csearch-min $1-csearch-min.mae
rm -f $1-csearch-min.com
rm -f $1-csearch-min.log
rm -f *mon.maegz
echo "DONE!"
date

# cleanup
cd ../../mae
rm -f $1.*
rm -f $1-mon.maegz
mv $1-* ../output/$1/
