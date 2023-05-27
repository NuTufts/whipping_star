#!/bin/bash

#SBATCH --job-name=sensband
#SBATCH --output=log-sensband-fixedue4
#SBATCH --partition batch
#SBATCH --time=1-00:00:00
#SBATCH --mem-per-cpu=10000
#SBATCH --array=0-999
#SBATCH --exclude=i2cmp001,i2cmp023,c1cmp003,c1cmp004,p1cmp075


CONTAINER=/cluster/tufts/wongjiradlabnu//larbys/larbys-container/singularity_minkowskiengine_u20.04.cu111.torch1.9.0_comput8.sif
TOP_DIR=/cluster/tufts/wongjiradlabnu/twongj01/oscfit/whipping_star/
SCRIPT_DIR=/cluster/tufts/wongjiradlabnu/twongj01/oscfit/whipping_star/gridscripts/
WORK_DIR=/cluster/tufts/wongjiradlabnu/twongj01/oscfit/whipping_star/out_dir/

mkdir -p $WORK_DIR

module load singularity/3.5.3
cvmfs_config probe fermilab.opensciencegrid.org uboone.opensciencegrid.org
singularity exec -B /cluster:/cluster,/cvmfs:/cvmfs ${CONTAINER} bash -c "cd ${SCRIPT_DIR} && source runjob_sens_band.sh ${TOP_DIR} ${WORK_DIR}"
