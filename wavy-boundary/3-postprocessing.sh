#!/bin/bash

echo "Activating environment..."
source postpro/bin/activate

echo "Performing postprocessing..."
python3 postprocessing.py
deactivate