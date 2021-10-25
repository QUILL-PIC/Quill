#!/bin/bash

# we do not want the script to continue if anything fails
set -e

YELLOW='\033[0;33m'
NOCOLOR='\033[0m'

config_filepath="$1"
if [ ! -f "$1" ]; then
    old_config_filepath=../quill3d-conf/quill.conf"$1"
    if [ -f $old_config_filepath ]; then
        echo -e "${YELLOW}Defaulting to config [$old_config_filepath]. Prefer using absolute paths instead! ${NOCOLOR}"
        config_filepath=$old_config_filepath
    else
        echo "Config file $1 not found!"
        exit 1
    fi
fi

# Create data folder if not exist; copy config file into it
folder=`./parse.py "$config_filepath" | grep -A2 data_folder | tail -n1`
if [[ -z "$folder" ]]
then
    folder="results"
fi
echo "Data folder: $folder"
sleep 1.5  # so the user has time to see the data folder

# Determine the number of threads specified in the config
threads=`./parse.py "$config_filepath" | grep -A1 n_sr | tail -n1`
if [[ -z "$threads" ]]
then
    threads=8
fi
echo "Threads: $threads"

mkdir -p $folder

cp "$config_filepath" $folder/

# 1. Parsing config file (parse.sh)
# 2. Running quill with input from step 1
# 3. Adding timestamps to the Quill's console output
# 4. Duplicating output to a log file (quill_log.txt)
mpirun -n $threads bash -c "stdbuf -o 0 ./parse.py \"$config_filepath\" | { ./quill; } 2>&1" | awk '{ print strftime("%Y-%m-%d %H:%M:%S\t"), $0; fflush(); }' | tee $folder/quill_log.txt
