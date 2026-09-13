#!/bin/bash

echo "Getting Basilisk..."
wget https://basilisk.fr/basilisk/basilisk.tar.gz
tar xzf basilisk.tar.gz

echo "Installing..."
cd basilisk/src
ln -s config.gcc config
make -k
make

echo "Adding qcc executable to PATH..."
echo "export BASILISK=$PWD" >> ~/.bashrc
echo 'export PATH=$PATH:$BASILISK' >> ~/.bashrc

echo "Installing additional packages..."
sudo apt install imagemagick ffmpeg

echo "Removing tar..."
rm basilisk.tar.gz