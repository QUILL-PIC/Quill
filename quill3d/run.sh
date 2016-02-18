#!/bin/bash

# Create data folder if not exist; copy config file into it
folder=`./parse.sh "$1" | grep -A2 data_folder | tail -n1`
echo "Data folder: $folder"
sleep 1.5  # so the user has time to see the data folder

mkdir -p $folder
cp ../quill3d-conf/quill.conf"$1" $folder/quill.conf"$1"

./parse.sh "$1" | ./quill
