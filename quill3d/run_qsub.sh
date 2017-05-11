#!/bin/sh
#PBS -N quill
#PBS -k oe
#PBS -l nodes=1:agolovanov:ppn=10
#PBS -l walltime=96:00:00
#PBS -q agolovanov

cd $PBS_O_WORKDIR
./run.sh $1
