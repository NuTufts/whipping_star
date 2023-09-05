#!/bin/bash

#SBATCH --job-name=fc2dprofile
#SBATCH --output=log-fc2dprofile-dm2-um4-set2
#SBATCH --partition batch,preempt,wongjiradlab
#SBATCH --time=8:00:00
#SBATCH --mem-per-cpu=8000
#SBATCH --array=0-999
##SBATCH --exclude=i2cmp001,i2cmp023,c1cmp003,c1cmp004


CONTAINER=/cluster/tufts/wongjiradlabnu//larbys/larbys-container/singularity_minkowskiengine_u20.04.cu111.torch1.9.0_comput8.sif
TOP_DIR=/cluster/tufts/wongjiradlabnu/twongj01/oscfit/whipping_star/
SCRIPT_DIR=/cluster/tufts/wongjiradlabnu/twongj01/oscfit/whipping_star/gridscripts/

# Ue4-Um4
#WORK_DIR=/cluster/tufts/wongjiradlabnu/twongj01/oscfit/whipping_star/workdir_mayagrid_fc2dprofile_ue4_um4_katiefitbounds/
#RUNLIST=${WORK_DIR}/rerun_w_starts.txt

# DM2-Um4
WORK_DIR=/cluster/tufts/wongjiradlabnu/twongj01/oscfit/whipping_star/workdir_mayagrid_fc2dprofile_dm2_um4_katiefitbounds/
RUNLIST=${WORK_DIR}/rerun_list_dm2_um4.txt

SCRIPTS=runjob_fc_2dmarginalize.sh
REDUCE_DIM=1 # remove Ue4
OFFSET=0

mkdir -p $WORK_DIR

module load singularity/3.5.3
singularity exec -B /cluster:/cluster,/cvmfs:/cvmfs ${CONTAINER} bash -c "cd ${SCRIPT_DIR} && source ${SCRIPTS} ${TOP_DIR} ${WORK_DIR} ${RUNLIST} ${REDUCE_DIM} ${OFFSET}"
