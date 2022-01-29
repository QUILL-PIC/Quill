---
layout: default
title: Run
nav_order: 2
---

# Input files

Example input files can be found in the `quill3d-conf/examples` folder.
All possible parameters in the input file are described in the `quill.conf.example` file.

# Run

To run Quill, use the `run.sh` script from the `quill3d` folder and pass the name of the input file as a parameter.
For example,
```
./run.sh /home/user/my-problem
```

By default, `run.sh` relies on `mpirun`.
Rewrite the script itself if different behavior is required.

The number of MPI threads is determined by the `n_sr` parameter in the input file and should never be explicitly set to a different value.

*Note:* In order to run Quill, a Python interpreter is required on all nodes;
it is used to parse the input file.