#!/bin/bash
#SBATCH --mail-user=dan.zhou@vumc.org
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --constraint=haswell|skylake #sandybridge|
#SBATCH --time=05:00:00
#SBATCH --mem-per-cpu=20G
#SBATCH --job-name=process_variants_and_annotate
#SBATCH --array=1-138
##SBATCH --account=xxxxxxx_lab

#REPEATED_JOBS=`cat /home/zhoud2/tmp/resubmit.list | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`    #only used to submit failed jobs

#ml GCC OpenMPI R
ml GCC/8.2.0  OpenMPI/3.1.4 Intel/2019.1.144  IntelMPI/2018.4.274 R/3.6.0

#ml GCC/5.4.0-2.26  OpenMPI/1.10.3 pandas/0.18.1-Python-2.7.12 numpy/1.11.1-Python-2.7.12 scipy/0.17.0-Python-2.7.12 R


Rscript /gpfs52/home/zhoud2/script/common_rare/common/step2_process_variants_and_annotate.r  ${SLURM_ARRAY_TASK_ID}

