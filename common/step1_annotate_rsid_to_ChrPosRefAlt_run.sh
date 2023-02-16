#!/bin/bash
#SBATCH --mail-user=dan.zhou@vumc.org
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --constraint=haswell|skylake #sandybridge|
#SBATCH --time=08:00:00
#SBATCH --mem-per-cpu=80G
#SBATCH --job-name=annotate
##SBATCH --account=xxxxxxx_lab
#SBATCH -o /data/xxxxxxx_lab/zhoud2/common_rare/common/magma_result/tmp_rownames/annotate_rsid_to_ChrPosRefAlt.log

#REPEATED_JOBS=`cat /home/zhoud2/tmp/resubmit.list | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`    #only used to submit failed jobs

ml GCC OpenMPI R
#ml GCC/8.2.0  OpenMPI/3.1.4 Intel/2019.1.144  IntelMPI/2018.4.274 R/3.6.0

#ml GCC/5.4.0-2.26  OpenMPI/1.10.3 pandas/0.18.1-Python-2.7.12 numpy/1.11.1-Python-2.7.12 scipy/0.17.0-Python-2.7.12 R

input=/data/xxxxxxx_lab/zhoud2/common_rare/common/magma_result/tmp/46_irnt.gwas.imputed_v3.both_sexes.tsv.gz
### output file name don't end with .gz, R script will gzip the output file
output=/data/xxxxxxx_lab/zhoud2/common_rare/common/magma_result/tmp_rownames/46_irnt.gwas.imputed_v3.both_sexes_add_rsid.txt

############ select the chr:pos:ref:alt colums index from input data, or one column contains chr:pos:ref:alt
############ if input_col_select only have one colunm, tell the sep="_" or ":", like chr1:2:A:t or 2:300:C:t
############ tell your input data is "b37" or "b38" version
Rscript /home/zhoud2/script/fun/annotate_rsid_to_ChrPosRefAlt.r -c 1 -s : -v b37 -i ${input} -o ${output}

