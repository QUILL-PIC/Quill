#!/bin/bash

cp ../quill3d-conf/quill.conf"$1" ./parse_temp

# Create data folder if not exist; copy config file into it
folder=` python ./parse.py | grep -A2 data_folder | tail -n1`
echo "Data folder: $folder"
sleep 1.5  # so the user has time to see the data folder

mkdir -p $folder
cp ../quill3d-conf/quill.conf"$1" $folder/quill.conf"$1"

# 1. Parsing config file (parse.sh)
# 2. Running quill with input from step 1
# 3. Adding timestamps to the Quill's console output
# 4. Duplicating output to a log file (quill_log.txt)
stdbuf -o 0 python ./parse.py | { ./quill; } 2>&1 | awk '{ print strftime("%Y-%m-%d %H:%M:%S\t"), $0; fflush(); }' | tee $folder/quill_log.txt

rm ./parse_temp
rm ./temp