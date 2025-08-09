# NekLB

The current version of NekLB is based on NekRS-24v

# Installation

Install the current version of NekRS:

```git clone https://github.com/Nek5000/NekRS.git nekrs_v24-dev```


```cd nekrs_v24-dev```

```git checkout v24-development```

Export the path:

```export NEKRS_HOME=<path-to-your-project>/.local/nekrs```

Installation:

```CC=mpicccc CXX=mpic++ FC=mpif77 ./build.sh -DCMAKE_INSTALL_PREFIX=$NEKRS_HOME```

   
# CPU:

compile the code by using Host:

```$NEKRS_HOME/bin/nrsmpi {name} {n} --backend serial --build-only 1```

# GPU:

compile the code by using Device:

```$NEKRS_HOME/bin/nrsmpi {name} {n} --backend CUDA --build-only 1```

# Run:

run the simulation by using n_processors, and the parameter in par file:

```mpirun -np {n_processors } nekrs --setup {.par}```

# Single phase
# Multi phase
