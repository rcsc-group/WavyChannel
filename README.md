## Running the simulation

The most recent version of the code is in `wavy-mpi.c`. To run, call `./mpirun-wavy.sh <npe>`, where `<npe>` is the number of desired MPI processes.

For the outputs, ensure you have the following folder structure to allow all outputs to be correctly saved (otherwise there will be an error):
```text
<repo>/
├── EBM_VOF/
├── output-mpi.h
├── wavy-mpi.c
├── mpirun-wavy.sh
└── fields/
    ├── interface/
    ├── level/
    ├── phase/
    └── pressure/
```

If there is an error pertaining to `virtual.h`, you may need to edit lines 26 and 27 of `grid/tree.h` from
```c
//#include "memindex/range.h"
#include "memindex/virtual.h"
```
to
```c
#include "memindex/range.h"
//#include "memindex/virtual.h"
```

## Problem setup

### Domain

 - Length: 1
 - Width: 0.01

### Boundary conditions

 - Left:
   - Prescribed pressure: $0.0$
   - Prescribed fluid: $1$
   - Free inflow
 - Right:
   - Prescribed pressure: $0.0$
   - Prescribed fluid: $0$
   - Free outflow
 - Top:
   - No slip
   - No flux
   - Contact angle: $60^\circ$
 - Bottom:
   - Symmetry
  
