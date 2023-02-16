
#ml GCC/8.2.0  OpenMPI/3.1.4 Intel/2019.1.144  IntelMPI/2018.4.274 R/3.6.0

library(data.table)
require(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)

info = as.data.frame(fread('/data/xxxxxxx_lab/zhoud2/common_rare/common/info/ukb_manually_merged.txt', header = T,stringsAsFactors = F))

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

#---------------

overlapped_phecode = intersect(common_trait_list, rare_trait_list)

ld_block = read.table('/data/coxvgi/zhoud2/anno/block/fourier_ls-all_hg19.bed', header = T,stringsAsFactors = F)
ld_block$chr = sapply(ld_block$chr, function(x) as.numeric(sub('chr','',x)))

ld_block_helper = function(chr, pos){
  ld_block_chr = ld_block[which(ld_block$chr == chr),]
  idx_1 = which(pos>ld_block_chr$start)
  idx_2 = which(pos<=ld_block_chr$stop)
  return(intersect(idx_1, idx_2))
}

#gencode
gencode = read.table('/data/coxvgi/zhoud2/anno/gencode/37/gencode.v32.GRCh37.txt',header = T,stringsAsFactors = F)
gencode = gencode[which(gencode$genetype %in% c('protein_coding','miRNA')),]
gencode$chr = as.numeric(sub('chr', '', gencode$chr))
gencode = gencode[!is.na(gencode$chr),]
gencode$pos = (gencode$left + gencode$right)/2
gencode = gencode[,c('genename','geneid','chr','pos')]
gencode$block = mapply(ld_block_helper, gencode[,'chr'], gencode[,'pos'])
gencode$block = paste0(gencode$chr, '_', gencode$block)
gencode = gencode[,c('genename','geneid','block')]



#---manhattan plot---
#---common---

manhattan_func = function(df, common_rare = 'common',tag = ' '){
  
  df$CHR = as.numeric(df$CHR)
  
  #pch
  df$pch = 'cir'
  top_gene_symbol = df_nondup[df_nondup[,paste0(common_rare,'_sig')] == 1,'gene_symbol']
  df[which(df$gene_symbol %in% top_gene_symbol),'pch'] = 'dim'
  df = df[order(df$pch, decreasing = T),]
  
  
  #label
  df$label = NA
  label_idx = which(df$gene_symbol %in% gene_list_label)
  df[label_idx, 'label'] = df[label_idx, 'gene_symbol']
  
  #horizontal line
  #logp_hline = -log(df[which.min(abs(df[,paste0('p_adj_',common_rare)]-0.05)),'P'],10)
  
  #pos
  df$plog = -log(df$P,10)
  df$pos=df$POS
  df$chr=df$CHR
  chr_box<-as.data.frame(matrix(data=NA,nrow=22,ncol=3))
  chr_df<-df[df$chr==1,]
  chr_box[1,1]<-max(chr_df$pos)
  chr_box[1,2]<-0
  chr_box[1,3]<-chr_box[1,1]/2
  
  for (i in 2:22){
    chr_df<-df[df$chr==i,]
    chr_box[i,1]<-max(chr_df$pos)
    chr_box[i,2]<-chr_box[i-1,1]+chr_box[i-1,2]
    df[which(df$chr==i),'pos']<-df[which(df$chr==i),'pos']+chr_box[i,2]
    chr_box[i,3]<-chr_box[i,2]+chr_box[i,1]/2
  }
  
  df$group<-ifelse(df$chr%%2==0,'odd','even')
  
  vline = df[label_idx,'pos']
  
  set.seed(1)
  manhattan <- ggplot(df, aes(x=pos, y=plog, label = label)) +
    scale_x_continuous(breaks=chr_box$V3, labels = c(as.character(seq(1,19)),' ','21',' ')) +  #x axis
    
    geom_vline(xintercept = vline, linetype="dashed", size = 0.1, color = 'grey')+
    #geom_hline(yintercept = logp_hline, linetype="dashed", size = 0.3, color = 'grey')+
    
    scale_y_continuous(trans='log10') +  #y axis log
    geom_point(data = df, aes(x = pos, y = plog, 
                              color = group,
                              alpha = pch,
                              shape = pch,
                              size = pch), show.legend = FALSE) + #points
    scale_color_manual(values = c("odd" = 'firebrick', "even" = 'dodgerblue4')) +
    scale_shape_manual(values = c("cir" = 20, 'dim' = 5))+
    scale_alpha_manual(values = c("cir" = 0.3, 'dim' = 0.95))+
    scale_size_manual(values = c("cir" = 0.75, 'dim' = 1.5))+
    
    labs(x = "Chromosome", y = "-log(P)") +
    theme(legend.position="none") +
    ggtitle(paste0(pheno_name,' (',common_rare,')'))+
    theme_bw() +  #rm background
    theme(panel.grid =element_blank())+ #rm grids
    
    labs(tag = tag)+
    
    
    geom_label_repel( #non overlapped labels
      size=2,
      #color = df$color,
      #nudge_x=5e4, #shift to the right
      segment.alpha = 0.2,  #transparent of segment
      min.segment.length = 1,
      segment.size = 0.2,
      fill=rgb(255, 255, 255, 210, maxColorValue=255),
      fontface='italic',
      seed = 2022
    )
  
  return(manhattan)
  
}



#generate figure

lay_out = function(...) {    
  x <- list(...)
  n <- max(sapply(x, function(x) max(x[[2]])))
  p <- max(sapply(x, function(x) max(x[[3]])))
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(n, p)))    
  
  for (i in seq_len(length(x))) {
    print(x[[i]][[1]], vp = grid::viewport(layout.pos.row = x[[i]][[2]], 
                                           layout.pos.col = x[[i]][[3]]))
  }
} 





#--qq plot--

dark_col = as.character(brewer.pal(12, "Paired")[c(2,8,4)])  #blue orange green

qqplot_func <- function(df) {
  #   
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  
  n = nrow(df)
  df_all = df
  df_05 = df[df$p_common<0.05,]
  df_01 = df[df$p_common<0.01,]
  df_001 = df[df$p_common<0.001,]
  
  #all
  
  tb <- data.frame(
    observed_all = -log10(sort(df_all$p_rare)),
    expected_all = -log10(ppoints(nrow(df_all))),
    
    clower   = -log10(qbeta(p = (1 - 0.95) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper   = -log10(qbeta(p = (1 + 0.95) / 2, shape1 = 1:n, shape2 = n:1))
  ) 
  
  
  tb$observed_05 = NA
  tb$expected_05 = NA
  tb[1:nrow(df_05),'observed_05'] = -log10(sort(df_05$p_rare))
  tb[1:nrow(df_05),'expected_05']  = -log10(ppoints(nrow(df_05)))
  
  tb$observed_01 = NA
  tb$expected_01 = NA
  tb[1:nrow(df_01),'observed_01'] = -log10(sort(df_01$p_rare))
  tb[1:nrow(df_01),'expected_01']  = -log10(ppoints(nrow(df_01)))
  
  tb$observed_001 = NA
  tb$expected_001 = NA
  tb[1:nrow(df_001),'observed_001'] = -log10(sort(df_001$p_rare))
  tb[1:nrow(df_001),'expected_001']  = -log10(ppoints(nrow(df_001)))
  
  qqplot = ggplot(tb) +
    
    #all
    geom_point(aes(expected_all, observed_all), shape = 1, size = 1.5) +
    
    #p<0.05
    geom_point(aes(expected_05, observed_05), shape = 1, size = 1.5, color = dark_col[1]) +
    
    #p<0.01
    geom_point(aes(expected_01, observed_01), shape = 1, size = 1.5, color = dark_col[2]) +
    
    #p<0.001
    geom_point(aes(expected_001, observed_001), shape = 1, size = 1.5, color = dark_col[3]) +
    
    annotate(geom="text", x=0.5, y=5, label="all")+
    annotate(geom="text", x=0.5, y=6, label=expression(paste("P"['common'], plain('<0.05'))), color = dark_col[1])+
    annotate(geom="text", x=0.5, y=7, label=expression(paste("P"['common'], plain('<0.01'))), color = dark_col[2])+
    annotate(geom="text", x=0.5, y=8, label=expression(paste("P"['common'], plain('<0.001'))), color = dark_col[3])+
    
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    geom_line(aes(expected_all, cupper), linetype = 2) +
    geom_line(aes(expected_all, clower), linetype = 2) +
    xlab(log10Pe) +
    ylab(log10Po) +
    ylim(0,10)+
    
    theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))
  
  return(qqplot)
}




#df

# result = as.data.frame(fread('/data/xxxxxxx_lab/zhoud2/common_rare/cor/result/Bonf.txt'))
# result = result[order(result$chi_test, decreasing = T), ]




phecode = '2296'	#Falls in the last year
phecode = '1239' #Current tobacco smoking
phecode = '189' #Townsend deprivation index at recruitment

phecode = '21001'  #BMI botton right
phecode = '20153'	 #Forced expiratory volume in 1-second (FEV1), predicted

phecode = '30080'  #Platelet count top right
phecode = '30790'  #Lipoprotein A top right
phecode = '30780'  #LDL top right
phecode = '30750'	#Glycated haemoglobin (HbA1c)
phecode = '30640' #Apolipoprotein B

#for(phecode in c('189','1239','2296','21001','20153','30080','30790','30780')){
#for(phecode in c('30750','30640')){

tag_i = 1

#for(phecode in c('30780','30750','21001','189')){
for(phecode in c('30690','30720','189','21001')){   #'Cholesterol' 'Cystatin C' 
  
  i = which(overlapped_phecode == phecode)  #30780LDL
  adj = 'FDR'
  evo = FALSE
    
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
  
  # p=rare$Pvalue_SKAT
  # z = qnorm(p/ 2)
  # lambda = round(median(z^2, na.rm = TRUE) / 0.454, 3)
  # 
  # p=rare$Pvalue_Burden
  # z = qnorm(p/ 2)
  # lambda = round(median(z^2, na.rm = TRUE) / 0.454, 3)

  rare$p = ifelse(rare$Pvalue_Burden > rare$Pvalue_SKAT, rare$Pvalue_SKAT, rare$Pvalue_Burden)
  rare = rare[order(rare$p,decreasing = T),]
  rare = rare[!duplicated(rare$gene_symbol),]
  rare = rare[!is.na(rare$p),]
  rare = rare[,c('gene_symbol','p','Pvalue_SKAT','Pvalue_Burden')]
  colnames(rare)[2] = 'p_rare'
  
  #merge
  df = merge(common, rare, by='gene_symbol')
  
  #rare variant-based clumping
  df = merge(df, gencode, by = 1)
  
  #df = df[order(df$p_rare),]
  df = df[order(df$p_rare*df$p_common),]
  
  #genes need a label (sig in both common and rare)
  df_nondup = df[!duplicated(df$block),]
  df_nondup$common_sig = ifelse(rank(df_nondup$p_common)<=100,1,0)
  df_nondup$rare_sig = ifelse(rank(df_nondup$p_rare)<=100,1,0)
  df_nondup$sig = df_nondup$common_sig * df_nondup$rare_sig
  
  common_sig_gene_symbol = df_nondup[df_nondup$common_sig == 1,'gene_symbol']
  rare_sig_gene_symbol = df_nondup[df_nondup$rare_sig == 1,'gene_symbol']
  
  gene_list_label = df_nondup[which(df_nondup$sig==1),'gene_symbol']
  
  
  #pheno info
  info = as.data.table(fread('/data/xxxxxxx_lab/zhoud2/common_rare/common/info/ukb_manually_merged.txt'))
  info = info[,c(1,2,3,4)]
  info$n_eff = ifelse(info$n_controls == 0, 
                      info$n_cases,
                      4*info$n_controls*info$n_cases / (info$n_controls+info$n_cases))
  
  pheno_name = info[which(info$phenocode == phecode), 'description']
  
  
  
  #
  df$pheno_name = as.character(pheno_name)
  write.table(df, paste0('/data/xxxxxxx_lab/zhoud2/common_rare/cor/eg/',phecode,'.txt'), quote = F, sep='\t',row.names = F)
  
  #QQ plot
  # qqplot = qqplot_func(df)
  # 
  # pdf(paste0('/data/xxxxxxx_lab/zhoud2/common_rare/cor/eg/QQ_',pheno_name,'.pdf'),width = 7, height = 7) #type = 'Xlib'
  # print(qqplot)
  # dev.off()
  
  
  
  df_common = df
  df_common$P = df$p_common
  df_common = df_common[df_common$P<0.1,]
  df_rare = df 
  df_rare$P = df$p_rare
  df_rare = df_rare[df_rare$P<0.1,]
  
  
  # common_fig = manhattan_func(df_common, common_rare = 'common',tag = letters[tag_i])
  # rare_fig = manhattan_func(df_rare, common_rare = 'rare', tag = ' ')
  # 
  # options(bitmapType='cairo')
  
  # pdf(paste0('/data/xxxxxxx_lab/zhoud2/common_rare/cor/eg/',pheno_name,'.pdf'),width = 8, height = 10) #type = 'Xlib'
  # lay_out(list(common_fig,1,1),
  #         list(rare_fig,2,1))
  # dev.off()
  
  tag_i = tag_i + 1

}
































