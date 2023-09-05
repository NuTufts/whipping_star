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

RUNLIST=$3

REDUCE_DIM=$4 # remove 0:dm2 1:Ue$ 2:Um4

OFFSET=$5 # offset to add to slurm task id

export PATH=${TOP_DIR}/build/Appearence/:${PATH}

let lineno="${SLURM_ARRAY_TASK_ID} + ${OFFSET} + 1"

let arrayid=`sed -n ${lineno}p ${RUNLIST} | awk '{ print $1 }'`
let startid=`sed -n ${lineno}p ${RUNLIST} | awk '{ print $2 }'`
echo "Array ID (corresponds to pseudo-experiment we fit): ${arrayid} start at experiment=${startid}"


#print out file number for log
pwd
cd ${TOP_DIR}
# setup container for ubdl
source misc/setenv_tufts_container.sh

#go to your working directory and make a temp directory for this job
cd ${WORK_DIR}
JOBDIR=job_fc2dprofile_${arrayid}
mkdir -p $JOBDIR
cd ${WORK_DIR}/$JOBDIR
# clear out a job dir if it has stuff in it
#rm ${WORK_DIR}/${JOBDIR}/*

echo "TOP_DIR: ${TOP_DIR}"
ls
pwd
echo "Check for DL3plus1_FCwregen_2d"
which DL3plus1_FCwregen_2d

# RESULTS FOLDER
mkdir -p ${WORK_DIR}/fit_results/
mkdir -p ${WORK_DIR}/pseudoexp_data/
mkdir -p ${WORK_DIR}/logs/

DL3plus1_FCwregen_2d $REDUCE_DIM $arrayid ${startid}  &> ${WORK_DIR}/${JOBDIR}/log_fc2dprofile_${arrayid}.txt

cp chis_seek* ${WORK_DIR}/fit_results/
cp pseudoexp_* ${WORK_DIR}/pseudoexp_data/
#cp log_* ${WORK_DIR}/logs/

# remove temporary job directory
cd ${WORK_DIR}

# rm -r ${JOBDIR}

echo Finished Job ${arrayid}
