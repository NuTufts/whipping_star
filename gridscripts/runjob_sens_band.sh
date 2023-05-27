#!/bin/bash

# This runs inside the container and builds the container
# We assume you did the first-time setup already
# //////////////////////////////////////////////////
# These jobs make the data to build
# a sensitivty band. This involves running
# pseudo-experiments from the CV model
# through `DL3plus1_data_wregen` which
# fits to the data and calculates the likelihood
# at each grid point
# /////////////////////////////////////////////////

# TOP_DIR is whipping_star folder 
TOP_DIR=$1
# WORK_DIR is where job output folders are written
WORK_DIR=$2

export PATH=${TOP_DIR}/build/Appearence/:${PATH}

# change  (+0) to i.e. 2000, if you have already run one set of jobs
let arrayid="$SLURM_ARRAY_TASK_ID"
echo $arrayid
#print out file number for log
pwd
cd ${TOP_DIR}
# setup container for ubdl
source misc/setenv_tufts_container.sh

#go to your working directory and make a temp directory for this job
cd ${WORK_DIR}
JOBDIR=job_sensband_${arrayid}
mkdir -p $JOBDIR
cd ${WORK_DIR}/$JOBDIR

echo "TOP_DIR: ${TOP_DIR}"
ls
pwd
echo "Check for DL3plus1_data_wregen"
which DL3plus1_data_wregen

# RESULTS FOLDER
mkdir -p ${WORK_DIR}/sensband_results/
mkdir -p ${WORK_DIR}/logs/

DL3plus1_data_wregen $arrayid  > ${WORK_DIR}/logs/log_sensband_${arrayid}.txt
# Rename output data
mv gridpts_data.txt gridpts_data_${arrayid}.txt
mv chis_data.txt chis_data_${arrayid}.txt

cp gridpts_data_${arrayid}.txt ${WORK_DIR}/sensband_results/
cp chis_data_${arrayid}.txt ${WORK_DIR}/sensband_results/

# remove temporary job directory
cd ${WORK_DIR}
echo Finished Job ${arrayid}
