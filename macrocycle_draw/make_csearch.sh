#!/bin/bash

# makes copy of the conformational search template
# with appropriate name $1, which should not have an extension

inputName=$1"-min.mae"
outputName="..\/output\/"$1"\/"$1"-csearch.maegz"
sed 1s/.*/"$inputName"/ csearch.com | sed 2s/.*/"$outputName"/ > ./mae/$1.com
