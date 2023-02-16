#ml GCC/8.2.0  OpenMPI/3.1.4 Intel/2019.1.144  IntelMPI/2018.4.274 R/3.6.0

args = as.numeric(commandArgs(TRUE))  #51

#i_start = (args-1)*10+1
#i_end = min(args*10,1043)

library(data.table)
library(stringr)
library(dplyr)
phenocode_i<-args

#for(phenocode_i in c(i_start:i_end))
#{
  
##############  read download files
data_file<-as.data.frame(fread('/data/xxxxxxx_lab/zhoud2/common_rare/common/info/ukb_manually_merged.txt'))
data_file$filename<-str_split(data_file$wget," ", simplify =T)[,4]


#load heritable traits
h2_info = as.data.frame(fread('/data/xxxxxxx_lab/zhoud2/common_rare/cor/result/top100.txt'))
data_file = data_file[which(data_file$phenocode %in% h2_info$phenocode),]
data_file = data_file[!duplicated(data_file$phenocode),]

################# download data
raw_download_cmd<-data_file$wget[phenocode_i]
setwd("/data/xxxxxxx_lab/zhoud2/common_rare/common/magma_result/tmp/")
system(raw_download_cmd)

################# get raw data filename
#raw_data<-list.files("/data/xxxxxxx_lab/zhoud2/common_rare/common/magma_result/tmp/")
############### get a phenocode raw data filename
raw_data<-data_file$filename[phenocode_i]
############# change raw data file .bgz format to .gz format
raw_data_rename<-str_replace_all(raw_data,".bgz",".gz")
system(paste("mv",raw_data,raw_data_rename))


############### get SNP and P information from raw data
source<-paste0("/data/xxxxxxx_lab/zhoud2/common_rare/common/magma_result/tmp/",raw_data_rename)
input_data<-fread(source)

########### delete raw data
system(paste("rm",raw_data_rename))  


############## read the whole variant from  46_irnt raw data
whole_variant<-fread("/data/xxxxxxx_lab/zhoud2/common_rare/common/magma_result/tmp_rownames/46_irnt.gwas.imputed_v3.both_sexes_whole_variant.txt.gz")


############## read annotate data
annotate_data<-fread("/data/xxxxxxx_lab/zhoud2/common_rare/common/magma_result/tmp_rownames/46_irnt.gwas.imputed_v3.both_sexes_variant_rsid.txt.gz")

############## check whether the variant column of input_data and annotate_data were the same
if(length(unique(input_data$variant==whole_variant$variant))==1&unique(input_data$variant==whole_variant$variant))
{
  Sys.time()
  ########### get MAGMA_input_data
  MAGMA_input_data<-merge(input_data,annotate_data,by.x="variant",by.y="variant",sort=F)
  MAGMA_input_data<-MAGMA_input_data[MAGMA_input_data$minor_AF>=0.05 ,c("rs_id_dbSNP151_GRCh38p7","pval")] #"rs_id_dbSNP151_GRCh38p7","pval"
  Sys.time()
  colnames(MAGMA_input_data)<-c("SNP","P")
  
  ############ calculate n_eff for MAGMA
  if(data_file$n_controls[phenocode_i]>0)
  {
    n_eff<-round(4*data_file$n_controls[phenocode_i]*data_file$n_cases[phenocode_i]/(data_file$n_controls[phenocode_i]+data_file$n_cases[phenocode_i]))
  }else
  {
    n_eff<-data_file$n_cases[phenocode_i]
  }
  
  ############### write MAGMA_input_data
  write.table(MAGMA_input_data,paste0("/data/xxxxxxx_lab/zhoud2/common_rare/common/magma_result/tmp_MAGMA_input_maf05/",data_file$phenocode[phenocode_i],"_MAGMA_input_data.txt"), quote =F,sep = '\t',row.names = F)
  
  ####################
  magma_cmd<-paste0("magma --bfile /data/coxvgi/zhoud2/tools/MAGMA/ref_geno/g1000_eur --pval ","/data/xxxxxxx_lab/zhoud2/common_rare/common/magma_result/tmp_MAGMA_input_maf05/",data_file$phenocode[phenocode_i],"_MAGMA_input_data.txt N=",n_eff," --gene-annot /data/coxvgi/zhoud2/tools/MAGMA/anno/37.genes.annot --out /data/xxxxxxx_lab/zhoud2/common_rare/common/magma_result_maf05/",data_file$phenocode[phenocode_i])
  system(magma_cmd)
  print(paste0("the phenocode was done: ",data_file$phenocode[phenocode_i]))
  
}else
{
  print(paste0("the phenocode whose variants does not match 46_irnt: ",data_file$phenocode[phenocode_i]))
}


#}

