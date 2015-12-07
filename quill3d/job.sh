#PBS -S /bin/bash
cd /home/eugn/2015.07.10-300TW-proposal/quill3d/
mpirun -hostfile /home/eugn/2015.07.10-300TW-proposal/quill3d/hfile -np 15 /home/eugn/2015.07.10-300TW-proposal/quill3d/a.out
