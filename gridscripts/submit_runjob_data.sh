#!/bin/bash

#SBATCH --job-name=data
#SBATCH --output=log-data
#SBATCH --partition wongjiradlab,preempt
#SBATCH --time=1-00:00:00
#SBATCH --mem-per-cpu=8000
#SBATCH --array=0

CONTAINER=/cluster/tufts/wongjiradlabnu/larbys/larbys-container/singularity_minkowskiengine_u20.04.cu111.torch1.9.0_comput8.sif
TOP_DIR=/cluster/tufts/wongjiradlabnu/twongj01/oscfit/whipping_star
SCRIPT_DIR=${TOP_DIR}/gridscripts/
WORK_DIR=${TOP_DIR}/workdir/

module load singularity/3.5.3
#cvmfs_config probe fermilab.opensciencegrid.org uboone.opensciencegrid.org
singularity exec -B /cluster:/cluster,/cvmfs:/cvmfs ${CONTAINER} bash -c "cd ${SCRIPT_DIR} && source runjob_data_wregen.sh ${TOP_DIR} ${WORK_DIR}"
