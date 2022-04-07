#!/bin/bash

# msim launch script
# (c) Alan Holt 2018

set -e

# Configure these parameters:
COUNT=1

msg="random drift"

DESC="TEST_2.5e-05_ppm-3e-06_init-2000_cap-10000"

PARAM_FILE="param-${DESC}.json"
DDIR="." # dir for results


OPTS[0]="-c ${COUNT}"
OPTS[1]="--desc=${DESC}"
OPTS[2]="-p ${DDIR}/${PARAM_FILE}"
OPTS[3]="-m ${msg}"

#docker run -u mito -v `pwd`/opt:/opt -d ${IMAGE}:${tag}  ${OPTS[*]}

./msim.py ${OPTS[*]}
