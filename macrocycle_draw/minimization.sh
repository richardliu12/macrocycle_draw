#!/bin/bash

# converts mol2 file $1
# $1 should not have a file extension

date
/n/sw/schrodinger/utilities/structconvert -imol2 ./$1.mol2 -omae ./mae/$1.mae
cat ./mae/$1.mae
inputName=$1".mae"
outputName=$1".maegz"
sed 1s/.*/"$inputName"/ minimize.com | sed 2s/.*/"$outputName"/ > ./mae/$1.com
cd mae
/n/sw/schrodinger/macromodel -LOCAL -WAIT $1.com
mv $1.maegz $1.gz
gunzip $1.gz
mv $1 $1-min.mae
rm -f $1.com
rm -f $1.log
rm -f $1-mon.maegz
/n/sw/schrodinger/utilities/structconvert -imae ./$1-min.mae -omol2 ./$1-min.mol2
date
