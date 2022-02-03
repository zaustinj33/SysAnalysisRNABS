#!/bin/bash

# Usage: bash 01_submit.sh <File with list of data> <Absolute path to directory>
# replace sh .... with execution bash file, eg 02_adapterTrim.sh

cat $1 | while read LINE
do
	sh ./02_adapterTrim.sh $LINE $2
done
