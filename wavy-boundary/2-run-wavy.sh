#!/bin/bash

echo "Compiling..."
qcc wavy.c -o wavy -O2 -L$BASILISK/gl -lglutils -lfb_tiny -lm

if [[ "$?" -ne 0 ]]; then
    echo "Compilation error"
    exit 1
fi

echo "Running..."
./wavy

open flow.mp4


