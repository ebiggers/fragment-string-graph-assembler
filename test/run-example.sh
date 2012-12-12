#!/usr/bin/env bash

make GENOME=E_coli.fa   \
     SAMPLE_SIZE=500000 \
     COVERAGE=5         \
     READ_LEN=100       \
     MIN_OVERLAP_LEN=25 \
     2>&1 | tee assemble.log
