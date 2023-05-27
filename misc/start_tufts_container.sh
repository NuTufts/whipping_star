#!/bin/sh

CONTAINER=/cluster/tufts/wongjiradlabnu//larbys/larbys-container/singularity_minkowskiengine_u20.04.cu111.torch1.9.0_comput8.sif

module load singularity/3.5.3
singularity shell -B/cluster:/cluster ${CONTAINER} bash
