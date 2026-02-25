#!/bin/bash

if [ -z "$1" ]; then
    echo "Usage: source thefile.sh <file_name>"
    return 1
fi


FILE_NAME="$1"
cp "/exp/e961/data/SBC-25-daqdata/$FILE_NAME" "/exp/e961/app/users/runze/data/"
tar -xvf "/exp/e961/app/users/runze/data/$FILE_NAME" "/exp/e961/app/users/runze/data/"
