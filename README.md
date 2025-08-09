# NekLB
NekLB
Installation
git clone https://github.com/Nek5000/NekRS.git nekrs_v24-dev
cd nekrs_v24-dev
git checkout v24-development
export NEKRS_HOME=<path-to-your-project>/.local/nekrs
 CC=mpicccc CXX=mpic++ FC=mpif77 ./build.sh \
   -DCMAKE_INSTALL_PREFIX=$NEKRS_HOME
CPU:
$NEKRS_HOME/bin/nrsmpi {name} {n} --backend serial --build-only 1
GPU:
$NEKRS_HOME/bin/nrsmpi {name} {n} --backend CUDA --build-only 1
Run:
mpirun -np {n_processors } nekrs --setup {.par}
Single phase
Multi phase
