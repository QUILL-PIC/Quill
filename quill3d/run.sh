#!/bin/bash

# we do not want the script to continue if anything fails
set -e

# Create data folder if not exist; copy config file into it
folder=`./parse.py "$1" | grep -A2 data_folder | tail -n1`
echo "Data folder: $folder"
sleep 1.5  # so the user has time to see the data folder

# Determine the number of threads specified in the config
threads=`./parse.py "$1" | grep -A1 n_sr | tail -n1`
if [[ -z "$threads" ]]
then
    threads=1
fi
echo "Threads: $threads"

mkdir -p $folder

cp ../quill3d-conf/quill.conf"$1" $folder/

# 1. Parsing config file (parse.sh)
# 2. Running quill with input from step 1
# 3. Adding timestamps to the Quill's console output
# 4. Duplicating output to a log file (quill_log.txt)
mpirun -n $threads bash -c "stdbuf -o 0 ./parse.py \"$1\" | { ./quill; } 2>&1" | awk '{ print strftime("%Y-%m-%d %H:%M:%S\t"), $0; fflush(); }' | tee $folder/quill_log.txt
