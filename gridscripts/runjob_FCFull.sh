#!/bin/bash

# This runs inside the container and builds the container
# We assume you did the first-time setup already

TOP_DIR=$1
WORK_DIR=$2

# change  (+0) to i.e. 2000, if you have already run one set of jobs
let arrayid="$SLURM_ARRAY_TASK_ID"
echo $arrayid
#print out file number for log
echo $FILEID
pwd
cd ${TOP_DIR}
# setup container for ubdl
source /usr/local/root/root-6.22.06/bin/thisroot.sh

#go to your working directory and make a temp directory for this job
cd ${WORK_DIR}
mkdir job_${arrayid}
cd job_${arrayid}

echo
ls
pwd
#copy your script to temp directory
# this assumes you already made your cxx code
cp ${WORK_DIR}/DL3plus1_FCwregen .
let startid="$arrayid+0*1900"
# runyour command

# make log directory
mkdir -p ${WORK_DIR}/logs

# make Results folder
mkdir -p ${WORK_DIR}/FCresults

./DL3plus1_FCwregen ${startid}  > ${WORK_DIR}/logs/${arrayid}_log.txt
mv chis*txt ${WORK_DIR}/FCresults/.

# move outputs to dedicated directory
#rm ${WORK_DIR}/textfiles2/chis_${arrayid}.txt 
#mv chis*txt  ${WORK_DIR}/FCresults/.
# remove temporary job directory
cd ${WORK_DIR}
rm -r job_${arrayid}
echo Finished
