---
layout: default
title: Run
nav_order: 2
---

# Input files

Example input files can be found in the [`conf/example`](https://github.com/QUILL-PIC/Quill/tree/master/conf/example) folder.
All possible parameters in the input file are described in the [`quill.conf.example`](https://github.com/QUILL-PIC/Quill/blob/master/conf/example/quill.conf.example) file.

# Run

To run Quill, use the [`run.sh`](https://github.com/QUILL-PIC/Quill/blob/master/run.sh) script from the root folder and pass the name of the input file as a parameter.
For example,
```
./run.sh /home/user/my-problem
```

By default, `run.sh` relies on `mpirun`.
Rewrite the script itself if different behavior is required.

The number of MPI threads is determined by the `n_sr` parameter in the input file and should never be explicitly set to a different value.

**Note:** In order to run Quill, a Python interpreter is required on all nodes;
it is used to parse the input file.

Parallelization in Quill is optimized for processess with adjacent ranks being on the same node/socket, so mapping by socket is advised.
Binding processes to cores is also advised to make sure the memory is always allocated.
For OpenMPI, the corresponding configuration is `mpirun --map-by socket --bind-to core`.
The configuration can be changed in the `run.sh` script.

## Running with slurm

Slurm is a popular workload manager for clusters.
See [Slurm documentation](https://slurm.schedmd.com/documentation.html) for comprehensive details.

The most common way to schedule a Quill job is through [`sbatch`](https://slurm.schedmd.com/sbatch.html).
An example shell-script [`run_slurm_example.sh`](https://github.com/QUILL-PIC/Quill/blob/master/run_slurm_example.sh) can be found in the root of the project.
As `run.sh`, it accepts the path to an inpit file
```
sbatch -N 1 -n 32 run_slurm_example.sh /home/user/my-problem
```

The most commonnly used arguments taken by `sbatch` are:
* `-N`, `--nodes=`: number of nodes (`-N 2` for 2 nodes).
* `-n`, `--ntasks=`: number of tasks (`-n 32` for 32 tasks)
* `-o`, `--output=`: file to redirect the process output to, e.g. `slurm-%x.%j.out` (refer to [Slurm documentation](https://slurm.schedmd.com/sbatch.html#SECTION_%3CB%3Efilename-pattern%3C/B%3E) for patterns)

Arguments can be explicitly passed to `sbatch` or written in the shell-script itself.

One can also use [`salloc`](https://slurm.schedmd.com/salloc.html) and run Quill using `run.sh` from the interactive Slurm session.

**Note**: the number of allocated slurm tasks (`ntasks`) does not affect the number of processees spawned by Quill (determined by the `n_sr` parameter, see above).
You have to make sure they are the same value.

## Running multiple jobs based on a template

Refer to the description in the [`run_slurm_example.py`](https://github.com/QUILL-PIC/Quill/blob/master/run_slurm_example.py) script for running several Slurm jobs based on a input template.