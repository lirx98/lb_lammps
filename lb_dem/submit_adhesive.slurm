#!/bin/bash --login

#SBATCH --job-name=lb_lammps
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00

# Replace [budget code] below with your project code
#SBATCH --account=e283-ruixin
#SBATCH --partition=standard
#SBATCH --qos=standard

#Setup the batch environment
# module load epcc-job-env

#srun launches the parallel program based on the SBATCH options
#export LD_LIBRARY_PATH=/home/e283/e283/lirx/lammps_install/lib64:$LD_LIBRARY_PATH

# 打印调试信息，确认路径生效
echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"
ldd ./dem_main | grep lammps
srun -distribution=block:block --hint=nomultithread ./dem_main
