#!/bin/bash

echo "Found $# arguments"
if [ "$#" -lt "1" ] ; then
    echo "Usage: paintsample1by1.sh <batch number>"
    echo "<batch_number>: batch number"
    exit 0
fi
cd $PBS_O_WORKDIR

batch="$1"
module load tools finestructure/4.1.1
module load tools parallel/20200922
module load anaconda3/4.4.0
module load intel/perflibs
module load gcc
module load R/3.6.1

cat batch_files/batch_commands$batch.txt | parallel