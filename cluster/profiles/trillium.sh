# Trillium GPU subcluster: 4 H100 per node, whole-node jobs, 24 h limit; submit from trig-login01.
FARGO_CLUSTER_NAME=trillium
FARGO_MODULES="StdEnv/2023 cuda/12.6 openmpi/4.1.5"
FARGO_CUDA_ARCH=sm_90
