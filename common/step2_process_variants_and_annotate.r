######################## process variants and annotate



library(data.table)
library(stringr)
########### need to run in first
annotation_data<-fread("/data/xxxxxxx_lab/zhoud2/common_rare/common/magma_result/tmp_rownames/46_irnt.gwas.imputed_v3.both_sexes_add_rsid.txt.gz",select=c(1,12)) ##### select variant and pvalue from raw data with annotation

write.table(annotation_data,"/data/xxxxxxx_lab/zhoud2/common_rare/common/magma_result/tmp_rownames/46_irnt.gwas.imputed_v3.both_sexes_variant_rsid.txt", quote =F,sep = '\t',row.names = F)

system("gzip /data/xxxxxxx_lab/zhoud2/common_rare/common/magma_result/tmp_rownames/46_irnt.gwas.imputed_v3.both_sexes_variant_rsid.txt")

system("rm /data/xxxxxxx_lab/zhoud2/common_rare/common/magma_result/tmp_rownames/46_irnt.gwas.imputed_v3.both_sexes_add_rsid.txt.gz")



########################################## save the whole variant from raw data, nrow(annotation_data) != nrow(input_data), because we delete the variants with two or more rsid
input_data<-fread("/data/xxxxxxx_lab/zhoud2/common_rare/common/magma_result/tmp/46_irnt.gwas.imputed_v3.both_sexes.tsv.gz",select=c(1))############## read raw data
annotate_data<-fread("/data/xxxxxxx_lab/zhoud2/common_rare/common/magma_result/tmp_rownames/46_irnt.gwas.imputed_v3.both_sexes_variant_rsid.txt.gz") ########## read annotation data

 ########### need to run in first
write.table(input_data,"/data/xxxxxxx_lab/zhoud2/common_rare/common/magma_result/tmp_rownames/46_irnt.gwas.imputed_v3.both_sexes_whole_variant.txt", quote =F,sep = '\t',row.names = F)

########### need to run in first
system("gzip /data/xxxxxxx_lab/zhoud2/common_rare/common/magma_result/tmp_rownames/46_irnt.gwas.imputed_v3.both_sexes_whole_variant.txt")

