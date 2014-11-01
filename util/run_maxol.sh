#!/bin/sh

export MAXOL_OUT_PATH="./out"

NT_END=10
NT_OUT=1

rm -f ${MAXOL_OUT_PATH}/*
./maxol ${NT_END} ${NT_OUT}