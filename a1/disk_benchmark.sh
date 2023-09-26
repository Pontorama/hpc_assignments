#!/bin/bash

dir_to_cp = /usr/include
ssd_dir = /run/mount/scratch/hpcuser058
hdd_dir = include_copy

start_time = "$(date)"
for i in {1..10}
do
  cp -r $dir_to_cp $hdd_dir
done
stop_time = "$(date)"

echo "$(($stop_time - $start_time))"


start_time = "$(date)"
for i in {1..10}
do
  cp -r $dir_to_cp $ssd_dir
done
stop_time = "$(date)"
