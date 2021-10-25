#!/bin/sh
#SBATCH -p kost_group
#SBATCH -o slurm-%x.%j.out
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 24:00:00
#SBATCH -J quill

./run.sh $1
