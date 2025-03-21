#!/bin/bash
#SBATCH --nodes=6
#SBATCH --ntasks-per-node=128
#SBATCH --ntasks=768
#SBATCH --nodelist=node02,node03,node05,node06,node07,node08
#SBATCH -e %j
#SBATCH -J 15_friction_10MPa_2e22_damper

## module load openmpi/gcc/64/latest 
module load openmpi/gcc/64/4.1.2
module load slurm
# ASP_R="/home/900358141/dealii-candi/LPO_hack/build/release/./aspect"
ASP_R="/home/900358141/dealii-candi/LPO_hack/build/./aspect"

mpirun $ASP_R continental_extension.prm

