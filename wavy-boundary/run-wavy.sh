#!/bin/bash

echo "Compiling..."
qcc wavy2.c -o wavy2 -O2 -L$BASILISK/gl -lglutils -lfb_tiny -lm

if [[ "$?" -ne 0 ]]; then
    echo "Compilation error"
    exit 1
fi

echo "Running..."
./wavy2

open flow.mp4


