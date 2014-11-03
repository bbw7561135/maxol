#!/bin/sh

export MAXOL_OUT_PATH="./binout"

NT_END=10
NT_OUT=1

MAXOL_TXT_OUT_PATH="./txtout"

rm -f ${MAXOL_OUT_PATH}/*
rm -f ${MAXOL_TXT_OUT_PATH}/*
./maxol ${NT_END} ${NT_OUT} 0.8

for f in `ls $MAXOL_OUT_PATH | grep ........[BE]`
do
../util/bin2txt $MAXOL_OUT_PATH/$f > ${MAXOL_TXT_OUT_PATH}/$f	
done
