#!/bin/bash

echo "Compiling..."
CC99='mpicc -std=c99' qcc -Wall -O2 -D_MPI=1 wavy.c -o wavy-mpi -lm -L$BASILISK/gl -lglutils -lfb_tiny

if [[ "$?" -ne 0 ]]; then
    echo "Compilation error"
    exit 1
fi

echo "Running..."
mpirun -np $1 ./wavy-mpi

open flow.mp4


