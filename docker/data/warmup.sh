#! /bin/bash -
CORES=16

printf 'Running selftest...\n'
python deconvolute.py selftest

printf 'Copying precompiled code for %s threads...\n' $CORES
seq $CORES | xargs -n1 -I{} -P$CORES --process-slot-var=SLOT \
    bash -c 'rsync -a ~/theano/slot_single ~/theano/slot_$SLOT'

printf 'Warmup done.\n'
