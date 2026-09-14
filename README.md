## Installation of Basilisk

This section is only required if you do not already have Basilisk installed on your machine.

Details of the Basilisk installation can be found on the [Basilisk website](https://basilisk.fr/src/INSTALL). The provided script `0-install-basilisk.sh` performs all the commands listed below for simplicity, and should work on a UNIX-like system. The only requirements are a C99-compatible compiler, and a version of `make` compatible with GNU make.

In the directory that you would like to install Basilisk, run the commands
```
wget https://basilisk.fr/basilisk/basilisk.tar.gz
tar xzf basilisk.tar.gz
```
To compile, run the commands
```
cd basilisk/src
ln -s config.gcc config
make -k
make
```
Now to add `qcc` to your PATH, run the commands
```
echo "export BASILISK=$PWD" >> ~/.bashrc
echo 'export PATH=$PATH:$BASILISK' >> ~/.bashrc
```
You may also require FFmpeg and ImageMagick. To install these, simply type
```
sudo apt install imagemagick ffmpeg
```

## Preparing directories

Basilisk should now be installed. To check this, open a new terminal and navigate to any directory. Run the command
```
qcc --version
```
If instead you get an error message, consult the installation instructions at the link above.

Now, run the script `1-create-directories.sh`. This will create the `interface` and `movies` directories for interfacial and field data respectively to be saved into, if they do not exist already. Furthermore, this sets up the virtual environment `postpro` for the Python postprocessing, and installs the required Python libraries. In particular, these are
 - NumPy,
 - SciPy,
 - Pandas,
 - Matplotlib.

## Running the simulation

The most recent version of the code is in `wavy.c`. To run, call `./2-mpirun-wavy.sh <npe>`, where `<npe>` is the number of desired MPI processes. For serial, call only `./2-run-wavy.c`.

<!-- For the outputs, ensure you have the following folder structure to allow all outputs to be correctly saved (otherwise there will be an error):
```text
<repo>/
├── EBM_VOF/
├── output-mpi.h
├── wavy.c
├── 2-mpirun-wavy.sh
├── 2-run-wavy.sh
├── interface/
└── movies/
``` -->

If there is an error pertaining to `virtual.h` on compilation, you may need to edit the Basilisk source code. Run the command
```
cd $BASILISK
```
to change directory to the source. Then, you will need to edit lines 26 and 27 of `grid/tree.h` from
```c
//#include "memindex/range.h"
#include "memindex/virtual.h"
```
to
```c
#include "memindex/range.h"
//#include "memindex/virtual.h"
```

## Postprocessing and reproducibility

The postprocessing can now be run with the script `3-postprocessing.sh`. This calculates the fitted growth rate and uncertainty in the fit of the growth rate, outputting these to the terminal. A plot showing interface position with time, the fitted position, and the scaled residuals for late times is also saved to the file `img-growth-rate.pdf`.

If additional precision is required in the fitted growth rate, this can be done by changing lines 47 and 48 of `postprocessing.py`.

If the Basilisk code is unable to run, precomputed data is also provided in the `interfaces` folder.