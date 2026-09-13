#!/bin/bash

echo "Creating directories for outputs..."
mkdir interface
mkdir movies

echo "Creating environment for postprocessing..."
python3 -m venv postpro
source postpro/bin/activate

echo "Installing packages..."
pip3 install numpy scipy pandas matplotlib jupyterlab ipykernel
deactivate
