
#ml GCC/8.2.0  OpenMPI/3.1.4 Intel/2019.1.144  IntelMPI/2018.4.274 R/3.6.0

library(data.table)
require(dplyr)
library(foreach)   
library(doParallel)   
library(ggplot2)
library(RColorBrewer)
library(ggrepel)


#-----load burden/skat-o results from rare variant----
phecode_helper_rare = function(x){
  x = sub('....$','',x)
  x = ifelse(substr(x,nchar(x),nchar(x)) == '_',
             sub('.$','',x),
             x)
  return(x)
}

rare_file_list = dir('/data/xxxxxxx_lab/zhoud2/common_rare/priority_needs_ukbb/ukbb_WES/')
rare_trait_list = sapply(rare_file_list, function(x) phecode_helper_rare(x))
tail(rare_trait_list)


#-----load MAGMA results from common variant----
phecode_helper_common = function(x){
  x = substr(x, nchar(x)-2, nchar(x))
  return(x)
}

common_file_list = dir('/data/xxxxxxx_lab/zhoud2/common_rare/common/magma_result_maf05/')

common_file_index = which(sapply(common_file_list, function(x) phecode_helper_common(x)) == 'out')
common_file_list = common_file_list[common_file_index]

common_trait_list = sapply(common_file_list, function(x) sub('.genes.out','',x))

common_annotation = read.table('/data/coxvgi/zhoud2/tools/MAGMA/anno/NCBI37.3.gene.loc',header = F,stringsAsFactors = F)
common_annotation = common_annotation[,c(1,6)]
colnames(common_annotation) = c('magma_id','gene_symbol')

#---merge common and rare asso results---
overlapped_phecode = intersect(common_trait_list, rare_trait_list)


#---load ld block info--------
ld_block = read.table('/data/coxvgi/zhoud2/anno/block/fourier_ls-all_hg19.bed', header = T,stringsAsFactors = F)
ld_block$chr = sapply(ld_block$chr, function(x) as.numeric(sub('chr','',x)))

ld_block_helper = function(chr, pos){
  ld_block_chr = ld_block[which(ld_block$chr == chr),]
  idx_1 = which(pos>ld_block_chr$start)
  idx_2 = which(pos<=ld_block_chr$stop)
  return(intersect(idx_1, idx_2))
}

#----load gencode gene annotation-----
gencode = read.table('/data/coxvgi/zhoud2/anno/gencode/37/gencode.v32.GRCh37.txt',header = T,stringsAsFactors = F)
gencode = gencode[which(gencode$genetype %in% c('protein_coding','miRNA')),]
gencode$chr = as.numeric(sub('chr', '', gencode$chr))
gencode = gencode[!is.na(gencode$chr),]
gencode$pos = (gencode$left + gencode$right)/2
gencode = gencode[,c('genename','geneid','chr','pos')]
gencode$block = mapply(ld_block_helper, gencode[,'chr'], gencode[,'pos'])
gencode$block = paste0(gencode$chr, '_', gencode$block)
gencode = gencode[,c('genename','geneid','block')]


#-----load phenotype info-------
info = as.data.frame(fread('/data/xxxxxxx_lab/zhoud2/common_rare/common/info/ukb_manually_merged.txt'))
info = info[,c(1,2,3,6)]
colnames(info)[4] = 'description'
info$n_eff = ifelse(info$n_controls == 0, 
                    info$n_cases,
                    4*info$n_controls*info$n_cases / (info$n_controls+info$n_cases))
info = info[!duplicated(info$description),]
info = info[-which(sapply(info$description, function(x) length(grep("predict",x)))==1),]

#-----load h2 for common variant from Ben Neale's h2 estimation using ldsc-----
common_h2<-data.frame(fread('/data/coxvgi/zhoud2/data/ukbb/bn/h2/ukb31063_h2_all.02Oct2019.tsv'))
common_h2 = common_h2[common_h2$sex == 'both_sexes',]
common_h2 = common_h2[,c('phenotype','h2_liability','h2_p')]
common_h2 = common_h2[order(common_h2$phenotype),]
common_h2 = common_h2[!duplicated(common_h2$phenotype),]
common_h2[which(common_h2$h2_p>0.05),'h2_liability'] = 0
common_h2[which(common_h2$h2_liability<0),'h2_liability'] = 0
common_h2$phenotype = sapply(common_h2$phenotype, function(x) sub('_irnt','',x))
colnames(common_h2)[1] = 'phenocode'
info = merge(info, common_h2, by = 'phenocode')
info = info[which(info$h2_liability>0),]
dim(info)

#--------evolution-------
# evo_df = read.table('/data/coxvgi/zhoud2/projects/gtex/info/EvoStats_final2.txt',header = T,stringsAsFactors = F)
# evo_df = evo_df[!is.na(evo_df$Mouse_dN.dS),]
# conserved = evo_df[evo_df$Mouse_dN.dS<=0.1,]


#----CORAC-kappa function-------

kappa_func <- function(matrix){
  c = nrow(matrix)
  n = sum(matrix)
  sum_po = 0
  sum_pe = 0
  for (i in 1:c){
    sum_po = sum_po + matrix[i,i]
    row = sum(matrix[i,])
    col = sum(matrix[,i])
    sum_pe = sum_pe + row * col
    po = sum_po/n
    pe = sum_pe/(n*n)
  }
  # print(po)
  # print(pe)
  kappa = (po - pe)/(1-pe)
  return(kappa)
}

# k = matrix(c(1000,60,60,5), nrow = 2)
# k = matrix(c(41,3,4,27), nrow = 2)
# kappa_func(k)


kappa_helper = function(df){
  
  ans_table = table(df$common_sig, df$rare_sig)
  
  result = list()
  
  if(length(ans_table)==1){
    result$table_c0_r0 = ans_table[1]
    result$table_c1_r0 = 0.1
    result$table_c0_r1 = 0.1
    result$table_c1_r1 = 0.1
    result$kappa = NA
    
  }else if(length(ans_table) == 2){
    result$table_c0_r0 = ans_table[1]
    result$table_c1_r0 = ans_table[2]
    result$table_c0_r1 = 0.1
    result$table_c1_r1 = 0.1
    result$kappa = NA
    
  }else if(length(ans_table) == 4){
    result$table_c0_r0 = ans_table[1]
    result$table_c1_r0 = ans_table[2]
    result$table_c0_r1 = ans_table[3]
    result$table_c1_r1 = ans_table[4]
    
    kappa_ans = kappa_func(as.matrix(ans_table))
    result$kappa_ans = kappa_ans
  }
  
  return(result)
  
}

#---for test---
# i = which(overlapped_phecode == '21001')  #BMI
# i = which(overlapped_phecode == '50')  #height
# i = which(overlapped_phecode == '30780')  #30780LDL

#-----------

cover_func = function(i, evo=FALSE, adj = 'Bonf', top_rank = FALSE, top_rank_n = 100, top_rank_n_common = 100, top_rank_n_rare = 100){
  
  #common
  common = as.data.frame(fread(paste0('/data/xxxxxxx_lab/zhoud2/common_rare/common/magma_result_maf05/', overlapped_phecode[i], '.genes.out')))
  common = merge(common, common_annotation, by = 1)
  common = common[!is.na(common$P),]
  
  common$POS = (common$START + common$STOP)/2
  common = common[,c('gene_symbol','CHR','POS','P')]
  colnames(common)[4] = 'p_common'
  
  #rare
  rare_file_i = rare_file_list[which(rare_trait_list == overlapped_phecode[i])]
  
  rare = as.data.frame(fread(paste0('/data/xxxxxxx_lab/zhoud2/common_rare/priority_needs_ukbb/ukbb_WES/',rare_file_i)))
  rare = rare[rare$annotation == 'pLoF|missense|LC',]
  rare$p = ifelse(rare$Pvalue_Burden > rare$Pvalue_SKAT, rare$Pvalue_SKAT, rare$Pvalue_Burden)
  rare = rare[order(rare$p,decreasing = T),]
  rare = rare[!duplicated(rare$gene_symbol),]
  rare = rare[,c('gene_symbol','p')]
  colnames(rare)[2] = 'p_rare'

  #merge
  df = merge(common, rare, by='gene_symbol')
  df = merge(df, gencode, by = 1)
  if(evo){
    df = df[which(df$gene_symbol %in% conserved$GeneName),]
  }
  
  #clumping
  #df = df[order(df$p_rare),]    #rare-based clumping
  df = df[order(df$p_rare*df$p_common),]     #common*rare-based clumping
  df = df[!duplicated(df$block),]
  
  #adj common
  if(adj == 'FDR'){
    df$p_adj_common = p.adjust(df$p_common, method = 'BH')
  }else{
    df$p_adj_common = p.adjust(df$p_common)
  }
  
  #adj rare
  if(adj == 'FDR'){
    df$p_adj_rare = p.adjust(df$p_rare, method = 'BH')
  }else{
    df$p_adj_rare = p.adjust(df$p_rare)
  }
  
  #define significance
  if(top_rank){
    df$common_sig = ifelse(rank(df$p_common)<=top_rank_n_common,1,0)
    df$rare_sig = ifelse(rank(df$p_rare)<=top_rank_n_rare,1,0)
  }else{
    df$common_sig = ifelse(df$p_adj_common<0.05,1,0)
    df$rare_sig = ifelse(df$p_adj_rare<0.05,1,0)
  }
  
  df$sig = df$common_sig * df$rare_sig
  
  #save
  write.table(df, paste0('/data/xxxxxxx_lab/zhoud2/common_rare/cor/for_each_trait_maf05/',overlapped_phecode[i],'.txt'),quote = F,sep='\t',row.names = F)
  
  #output
  result = list()
  
  result$phecode = overlapped_phecode[i]
  
  #skip traits with less than 100 rare-based sig genes after FDR correction
  p_adj_rare_top_100 = df[100,'p_adj_rare']
  result_pruned = kappa_helper(df)
  result = c(result, result_pruned, p_adj_rare_top_100)
  
  return(result)
  
}

#-------run--------

registerDoParallel(cores = max(detectCores()-1, 1))

trials <- length(overlapped_phecode)

output <- foreach(i = 1:trials, .combine=rbind) %dopar% {
  ans = cover_func(i, adj = 'FDR', evo = FALSE, top_rank = TRUE)
  cat(paste0('INFO ',i,' completed \n'))
  unlist(ans)
}

stopImplicitCluster()

#for deBUG
# for(i in 1:length(overlapped_phecode)){
#   print(i)
#   tmp = chisq_func(i, adj = 'Bonf', evo = FALSE)
# }

output = as.data.frame(output,stringsAsFactors = F)
#output[,1] = as.character(output[,1])
colnames(output)[ncol(output)] = 'p_adj_rare_top_100'

for(j in seq(2,length(output))){
  output[,j] = as.numeric(output[,j])
}
df = merge(as.data.frame(info), output, by = 1)

# df$or = df$table_c1_r1 * df$table_c0_r0 / df$table_c1_r0 / df$table_c0_r1

write.table(df,'/data/xxxxxxx_lab/zhoud2/common_rare/cor/result/top100_maf05_kappa.txt',quote = F,sep='\t',row.names = F)



