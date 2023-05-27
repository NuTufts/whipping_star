#!/bin/bash

# This runs inside the container and builds the container
# We assume you did the first-time setup already
# //////////////////////////////////////////////////
# This script runs the code to make a delta-Chi2 surface
# on the real data.
#
# This involves running `DL3plus1_data_wregen` which
# calculates the best fit to the data along
# with the Chi2 at each grid point
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
JOBDIR=job_datafit_${arrayid}
mkdir -p $JOBDIR
cd ${WORK_DIR}/$JOBDIR
rm ${WORK_DIR}/${JOBDIR}/*.txt

echo "TOP_DIR: ${TOP_DIR}"
ls
pwd
echo "Check for DL3plus1_data_wregen"
which DL3plus1_data_wregen

# RESULTS FOLDER
DIR_RESULTS=${WORK_DIR}/results_datafit/
DIR_LOGS=${WORK_DIR}/logs_datafit/
mkdir -p ${DIR_RESULTS}
mkdir -p ${DIR_LOGS}

# Run with no argument. Will run grid over the data.
DL3plus1_data_wregen  > ${DIR_LOGS}/log_datafit_${arrayid}.txt

# Rename output data
mv gridpts_data.txt gridpts_data_${arrayid}.txt
mv chis_data.txt chis_data_${arrayid}.txt

cp gridpts_data_${arrayid}.txt ${DIR_RESULTS}/
cp chis_data_${arrayid}.txt ${DIR_RESULTS}/

# remove temporary job directory
cd ${WORK_DIR}
echo Finished Job ${arrayid}
