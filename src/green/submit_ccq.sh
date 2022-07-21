#!/bin/bash
## The computing resources we need
#SBATCH --nodes=4
##SBATCH --ntasks-per-node=28     ## For the [ccq], [gen] and [preempt] subcluster 
##SBATCH --ntasks-per-node=44     ## For the [ib] subcluster
##SBATCH --ntasks-per-node=40     ## For the [bnl] subcluster
##SBATCH --cpus-per-task=1 #for openmp

## The submitted partition of the job
#SBATCH --partition=ccq

#SBATCH -C opa

## Run this job on specific nodes
##SBATCH --nodelist=worker1035,worker1036

## Exclude specific nodes for this job
##SBATCH --exclude=worker1035,worker1036

## The limited computation time of the job
#SBATCH --time=7-00:00:00

## The memory usage of the job
##SBATCH --mem-per-cpu=1000mb

## The name of present job
#SBATCH --job-name=cs

## Send me a email if the job is finished
#SBATCH --mail-type=END
#SBATCH --mail-user=chenkun0228@gmail.com

## The standard output for the job
#SBATCH --output=slurm-%j.out

## The Error output file for the job
#SBATCH --error=slurm-%j.err

#mkdir data
#cp -r ./groups* ./data
#cp -r ./*.exe ./data
#cp parameter ./data
#cd data

## Get the information of tasks and nodes of present job
echo "#########################################################" >  host.txt
echo "SLURM_JOB_NUM_NODES  =" $SLURM_JOB_NUM_NODES               >> host.txt
echo "SLURM_JOB_NODELIST   =" $SLURM_JOB_NODELIST                >> host.txt
echo "SLURM_NTASKS         =" $SLURM_NTASKS                      >> host.txt
echo "SLURM_TASKS_PER_NODE =" $SLURM_TASKS_PER_NODE              >> host.txt
echo "#########################################################" >> host.txt

module load slurm
module load julia
module load openmpi4

## Just run the job
cd $SLURM_SUBMIT_DIR
#echo $SLURM_SUBMIT_DIR
/mnt/home/kunchen/.julia/bin/mpiexecjl julia /mnt/home/kunchen/project/EFT_UEG/green/greenMC.jl >> output.dat
#mpirun julia /mnt/home/kunchen/ceph/randomcircuit6/ConservedCircuit/evolve_2sector.jl
#mpirun julia evolve_ancilla.jl
#mkdir data
