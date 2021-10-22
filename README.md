# QUILL: 3D QED particle-in-cell code

**QUILL** (simulator for **QU**antum effects in **I**ntense **L**aser-p**L**asma  interactions) is a fully three-dimensional parallel particle-in-cell (PIC) code developed at the Institute of Applied Physics RAS, Nizhny Novgorod, Russia.
To our knowledge, it was the first PIC code with implementation of the Monte Carlo QED approach to investigate the development of electron–positron cascades.

The code is able to model the following processes using the Monte Carlo technique:

* photon emission by an electron in the strong field, with radiation reaction effects;
* electron–positron pair creation from gamma photons (Breit–Wheeler process).

The Maxwell solvers implemented in the code are FDTD, NDFX (the scheme used in A. Pukhov's VLPL code), and hybrid five-point FDTD (the scheme reduces numerical Cherenkov instability).
The particles pushers implemented in the code use Vay or Boris scheme.

# Dependencies

In order to build and run Quill, the following dependencies are required:
* C++ compiler with C++11 support;
* Make;
* MPI implementation;
* Python3 interpreter (*optional*: NumPy and Matplotlib for data analysis).

Quill is developed and tested on Linux with the use of *g++* and *clang++* compilers and *OpenMPI*.

# Build

To build Quill, run `make all` in the `quill3d` folder.

*Note*: Running `make` without target `all` does not build the required `chameleon` package.
Use `make` without targets only if `chameleon` is already built.

By default, the build process invokes `mpicxx` which should automatically provide necessary include and library paths for the compiler and linker.
Use the `CXX` environment variable to override this behavior.

When OpenMPI is used, the `OMPI_CXX` environment variable can be used to select the C++ compiler, e.g.
```
OMPI_CXX=clang++ make all
```

# Input files

Input files for Quill should be placed in the `quill3d-conf` folder.
An input file must begin with the `quill.conf` prefix.

Example input files can be found in the `quill3d-conf/examples` folder.
All possible parameters in the input file are described in the `quill.conf.example` file.

# Run

To run Quill, use the `run.sh` script from the `quill3d` folder and pass the name of the input file without the `quill.conf` prefix as a parameter.
For example,
```
./run.sh .my-problem
```
will run Quill using the `quill-conf/quill.conf.my-problem` input file.

By default, `run.sh` relies on `mpirun`.
Rewrite the script itself if different behavior is required.

The number of MPI threads is determined by the `n_sr` parameter in the input file and should never be explicitly set to a different value.

*Note:* In order to run Quill, a Python interpreter is required on all nodes;
it is used to parse the input file.

# Analyze results

Results can be graphically analyzed with the `qplot` Python package located in the `quill3d/python` folder.
Refer to the documentation within the package itself.
Documentation in IPython or Jupyter Notebooks/Lab can also be accessed using the `?` command, e.g. `qplot.density?`

Low-level data reading is available in the `resread` package.

The `qplot` package depends on the `matplotlib` and `numpy` packages; `resread` depends only on `numpy`.
