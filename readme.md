## do some physics note


## Quantum Espresso installation

[ref](https://pranabdas.github.io/espresso/setup/install/)

```linux
sudo apt update && sudo apt upgrade

sudo apt install --no-install-recommends libfftw3-dev quantum-espresso
```

## Quantum Espresso installation with sources 

[ref](https://pranabdas.github.io/espresso/setup/install/)

``` 
sudo apt install --no-install-recommends autoconf build-essential ca-certificates gfortran libblas3 libc6 libfftw3-dev libgcc-s1 liblapack-dev wget

sudo apt install --no-install-recommends libopenmpi-dev libscalapack-openmpi-dev libelpa17  #  use libelpa4 on Ubuntu 20.04

wget https://gitlab.com/QEF/q-e/-/archive/qe-7.2/q-e-qe-7.2.tar.gz

tar -zxvf q-e-qe-7.2.tar.gz

cd q-e-qe-7.2
./configure

# compile individual packages
make pw
# or compile everything
make all
# we can parallelize e.g., below command uses 4 CPUs
make -j4 all

# use the correct path if it differs from mine
echo 'export PATH="/root/q-e-qe-7.2/bin:$PATH"' >> ~/.bashrc

source ~/.bashrc

pdflatex filename.tex

pw.x -inp inputfile > outputfile
# For parallel version
mpirun -np 12 pw.x -inp inputfile > outputfile

make run-tests

sudo apt install tcl tcllib

wget "http://pwtk.ijs.si/download/pwtk-2.0.tar.gz"

tar -zxvf pwtk-2.0.tar.gz

echo 'export PATH="/root/pwtk-2.0:$PATH"' >> ~/.bashrc
source ~/.bashrc

```


