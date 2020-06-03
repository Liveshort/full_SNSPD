#!/bin/bash

# build
cd C/
make all

MKVAL="$?"
if [ "$MKVAL" -ne "0" ]
    then
        echo "Something went wrong compiling the program (error $MKVAL)... Exiting..."
        exit "$MKVAL"
fi

# determine number of physical cores, and put it in the right environment variable
if [ -z "$3" ]
    then
        echo "No number of threads given... Using nproc information..."
        NUM_PROCS=`nproc`
    else
        NUM_PROCS=$3
fi
export OMP_NUM_THREADS=$NUM_PROCS

# run
cd ../
time C/bld/simulation $1 $2

SIMVAL="$?"
if [ "$SIMVAL" -ne "0" ]
    then
        echo "Something went wrong running the program (error $SIMVAL)... Exiting..."
        exit "$SIMVAL"
fi

# plot
cd Python/
#python3 plot.py
cd ../
