
#ml GCC/8.2.0  OpenMPI/3.1.4 Intel/2019.1.144  IntelMPI/2018.4.274 R/3.6.0; R

#-----settings-----

args=as.numeric(commandArgs(TRUE))  #subjob id

run_i=1
run_array=list()
for (h2 in c(0.1,0.3,0.5,0.8)[2]){ #h2
  for (n in c(10000, 5000, 3000, 1000)){ 
    for (prop_causal_genes in c(0.1)){ 
      for (shuffle_rare in c(1, 2)){ 
        for (genetic_model in c(1,2,3,4)){ #genetic architecture
          for (seed_i in c(1:100)){  
            for (common_rare_sharing_proportion in c(0,0.5,1)[2]){
              run_array[[run_i]]=c(h2,n,prop_causal_genes,
                                   shuffle_rare,genetic_model,seed_i,
                                   common_rare_sharing_proportion)
              run_i=run_i+1
            }
          }
        }
      }
    }
  }
}

run_array = run_array[[args]]
h2 = run_array[1]
n = run_array[2]
prop_causal_genes = run_array[3]
shuffle_rare = switch(run_array[4], TRUE, FALSE)
genetic_model = switch(run_array[5],'a','b','c','d')
# a=Infinitesimal architecture
# b=with negative selection
seed_i = run_array[6]
pi0 = 0
common_rare_sharing_proportion = run_array[7]


library(data.table)


#-----load annotation-----
gencode = read.table('/data/***/anno/gencode/37/gencode.v32.GRCh37.txt',header = T,stringsAsFactors = F)
gencode = gencode[which(gencode$genetype %in% c('protein_coding')),]
gencode = gencode[gencode[,'chr'] == 'ukb', ]
gencode = gencode[,c('genename','geneid','left','right')]

#-----load snp info for whole chr 22-----
whole_chr_bim = as.data.frame(fread('/data/***/common_rare/simulation/ukb/ukb_whole/ukb_10k_maf0.001.bim'))
whole_chr_bim[,'pos_id'] = as.numeric(row.names(whole_chr_bim))
snp_with_rsid_index = which(substr(whole_chr_bim[,'V2'], 1, 2) == 'rs')
whole_chr_bim = whole_chr_bim[snp_with_rsid_index,]


#-----load dosage for whole chr 22-----
whole_chr_dosage = as.data.frame(fread('/data/***/common_rare/simulation/ukb/ukb_whole/ukb_10k_maf0.001.raw'))
id_df = whole_chr_dosage[,c('FID','IID')]
whole_chr_dosage = whole_chr_dosage[,c(whole_chr_bim[,'pos_id']+6)]
colnames(whole_chr_dosage) = whole_chr_bim$V2
whole_chr_bim[,'pos_id'] = seq(1, nrow(whole_chr_bim)) #reset


#-----define a set of causal genes-----
gene_list = dir('/data/***/common_rare/simulation/ukb/ukb_dosage_each_gene/')
n_causal_genes = round(length(gene_list)*prop_causal_genes)

#common causal gene
set.seed(seed_i)
common_causal_gene = sample(gene_list, n_causal_genes, replace = F)

#rare causal gene
n_shared_genes = round(common_rare_sharing_proportion * n_causal_genes)
rare_causal_gene = sample(common_causal_gene, n_shared_genes, replace = F)
n_rare_causal_gene_append = n_causal_genes - length(rare_causal_gene)
rare_causal_gene_append = sample(setdiff(gene_list, common_causal_gene), n_rare_causal_gene_append, replace = F)
rare_causal_gene = c(rare_causal_gene, rare_causal_gene_append)
# length(intersect(rare_causal_gene, common_causal_gene))


#-----load the snps-----

#load maf
freq = as.data.frame(fread('/data/***/common_rare/simulation/ukb/ukb_whole/ukb_10k_maf0.001.frq'))
freq = merge(freq, whole_chr_bim[,c('V2','V4','pos_id')], by.x = 'SNP', by.y = 'V2')
colnames(freq) = c("SNP", "CHR", "A1", "A2", "MAF", "NCHROBS", "BP", "pos_id")
common_variant_df = freq[freq[,'MAF']>=0.01, ]
rare_variant_df = freq[freq[,'MAF']<0.01, ]


#-----shuffle rare-----
if(shuffle_rare){
  whole_chr_dosage = readRDS('/data/***/common_rare/simulation/ukb/ukb_dosage_each_gene_rare_shuffled_all/ukb_dosage_each_gene_rare_shuffled.rds')
}



#func for extract_snps

causal_gene_list = common_causal_gene
variant_df = common_variant_df

extract_snps = function(causal_gene_list, variant_df){
  
  gencode_tmp = gencode[which(gencode[,'geneid'] %in% causal_gene_list),]
  
  i = 1
  for(i in 1:nrow(gencode_tmp)){
    geneid = gencode_tmp[i,'geneid']
    left = gencode_tmp[i,'left']
    right = gencode_tmp[i,'right']
    
    index = intersect(which(variant_df[,'BP']>=left), 
                      which(variant_df[,'BP']<=right))
    snp_tmp = variant_df[index, ]
    
    if(i == 1){
      snp_df = snp_tmp
    }else{
      snp_df = rbind(snp_df, snp_tmp)
    }
  }
  return(snp_df)
}


snp_df_common = extract_snps(causal_gene_list = common_causal_gene,
                             variant_df = common_variant_df)
snp_df_common = snp_df_common[!duplicated(snp_df_common[,'SNP']),]
snp_df_rare = extract_snps(causal_gene_list = rare_causal_gene,
                           variant_df = rare_variant_df)
snp_df_rare = snp_df_rare[!duplicated(snp_df_rare[,'SNP']),]
snp_df = rbind(snp_df_common, snp_df_rare)
snp_df = snp_df[!duplicated(snp_df[,'SNP']), ]



#-----effect size-----

#load pruned snps
pruned_bim = as.data.frame(fread('/data/***/common_rare/simulation/ukb/ukb_whole/ukb_10k_maf0.001_pruned.bim'))

snp_pruned_df = snp_df[which(snp_df[,'SNP'] %in% pruned_bim[,'V2']),]

i = 1

N_causal_SNPs = nrow(snp_pruned_df)
set.seed(seed_i+10)

if(genetic_model == 'a'){
  alpha = 0
  
}else if(genetic_model == 'b'){
  alpha = -0.37
  
}else if(genetic_model == 'c'){
  alpha = -1
  
}else if(genetic_model == 'd'){
  alpha = 0.5
}

snp_pruned_df[,'effect'] = sapply(snp_pruned_df$MAF, function(f) rnorm(1,mean=0,sd=((f*(1-f))^(1+alpha))^0.5))
k_constant=(h2/N_causal_SNPs/var(snp_pruned_df[,'effect']))
snp_pruned_df[,'effect'] = snp_pruned_df[,'effect']*k_constant^0.5*(1-pi0)



#-----define phenotype-----

df = whole_chr_dosage[,snp_pruned_df[,'pos_id']]
df = apply(df, MARGIN = 2, function(x) ifelse(is.na(x), median(x, na.rm = T), x)) #impute
df = apply(df, MARGIN = 2, function(x) scale(x)) #scale

common_index = which(colnames(df) %in% common_variant_df[,'SNP'])
rare_index = which(colnames(df) %in% rare_variant_df[,'SNP'])

df_common = df[, common_index]
df_rare = df[, rare_index]

snp_df_common = snp_pruned_df[sapply(colnames(df_common), function(x) which(snp_pruned_df[,'SNP'] == x)[1]),]
snp_df_rare = snp_pruned_df[sapply(colnames(df_rare), function(x) which(snp_pruned_df[,'SNP'] == x)[1]),]

# table(colnames(df_common) == snp_df_common[,'SNP'])

score_common = as.numeric(as.matrix(df_common) %*% as.matrix(snp_df_common[,'effect']))
score_rare = as.numeric(as.matrix(df_rare) %*% as.matrix(snp_df_rare[,'effect']))
score = score_common+score_rare

set.seed(i+2023)
score_scaled = scale(score)*h2 + rnorm(length(score), mean = 0, sd = (1-h2)^0.5)


#-----sample index-----
set.seed(i+100)
sample_index = sample(nrow(whole_chr_dosage), n, replace = F)

#-----run gwas for common variants-----
trait_df = cbind(id_df, as.data.frame(score_scaled))
colnames(trait_df)[3] = 'trait'
write.table(trait_df[sample_index,], paste0("/data/***/common_rare/simulation/common_magma/trait/h2_",h2,"_geneticmodel_",genetic_model,"_shareprop_",common_rare_sharing_proportion,"_shufflerare_",shuffle_rare,"_n_",n,"_propcausal_",prop_causal_genes,"_seed_",seed_i,".txt"), quote = F, row.names = F, sep='\t')

cmd = paste0("plink --bfile /data/***/common_rare/simulation/ukb/ukb_whole/ukb_10k_maf0.001 --pheno  /data/***/common_rare/simulation/common_magma/trait/h2_",h2,"_geneticmodel_",genetic_model,"_shareprop_",common_rare_sharing_proportion,"_shufflerare_",shuffle_rare,"_n_",n,"_propcausal_",prop_causal_genes,"_seed_",seed_i,".txt --pheno-name trait --allow-no-sex --maf 0.01 --linear --out /data/***/common_rare/simulation/common_magma/asso_plink/h2_",h2,"_geneticmodel_",genetic_model,"_shareprop_",common_rare_sharing_proportion,"_shufflerare_",shuffle_rare,"_n_",n,"_propcausal_",prop_causal_genes,"_seed_",seed_i)
system(cmd, wait = T)

#-----MAGMA test for each gene-----
magma_cmd<-paste0("magma --bfile /data/***/common_rare/simulation/ukb/ukb_whole/ukb_10k_maf0.001 --pval ","/data/***/common_rare/simulation/common_magma/asso_plink/h2_",h2,"_geneticmodel_",genetic_model,"_shareprop_",common_rare_sharing_proportion,"_shufflerare_",shuffle_rare,"_n_",n,"_propcausal_",prop_causal_genes,"_seed_",seed_i,".assoc.linear N=",length(sample_index)," --gene-annot /data/***/tools/MAGMA/anno/37.genes.annot --out /data/***/common_rare/simulation/common_magma/asso_plink/h2_",h2,"_geneticmodel_",genetic_model,"_shareprop_",common_rare_sharing_proportion,"_shufflerare_",shuffle_rare,"_n_",n,"_propcausal_",prop_causal_genes,"_seed_",seed_i,".txt")
system(magma_cmd)



#-----SKAT test for each gene-----

library(SKAT)

skat_func = function(geneid){
  
  if(shuffle_rare){
    dosage_rare = as.data.frame(fread(paste0('/data/***/common_rare/simulation/ukb/ukb_dosage_each_gene_rare_shuffled/',geneid)))
  }else{
    dosage_rare = as.data.frame(fread(paste0('/data/***/common_rare/simulation/ukb/ukb_dosage_each_gene/',geneid)))
  }
  
  rare_index = which(colnames(dosage_rare) %in% rare_variant_df[,'SNP'])
  dosage_rare = dosage_rare[,rare_index]
  
  obj = SKAT_Null_Model(score[sample_index] ~ 1, out_type="C")
  SKAT_P = SKAT(as.matrix(dosage_rare[sample_index,]), obj)$p.value
  
}

skat_rs = sapply(gene_list, function(x) skat_func(x))

skat_rs = data.frame(geneid = gene_list,
                     p_common = skat_rs)

#-----load results-----

magma_rs = as.data.frame(fread(paste0("/data/***/common_rare/simulation/common_magma/asso_plink/h2_",h2,"_geneticmodel_",genetic_model,"_shareprop_",common_rare_sharing_proportion,"_shufflerare_",shuffle_rare,"_n_",n,"_propcausal_",prop_causal_genes,"_seed_",seed_i,".txt.genes.out")))

common_annotation = read.table('/data/***/tools/MAGMA/anno/NCBI37.3.gene.loc',header = F,stringsAsFactors = F)
common_annotation = common_annotation[,c(1,6)]
colnames(common_annotation) = c('magma_id','gene_symbol')

magma_rs = merge(magma_rs, common_annotation, by = 1)

gencode = read.table('/data/***/anno/gencode/37/gencode.v32.GRCh37.txt',header = T,stringsAsFactors = F)
gencode = gencode[,c('geneid','genename')]
magma_rs = merge(magma_rs, gencode, by.x = 'gene_symbol', by.y = 'genename')


rs = merge(magma_rs, skat_rs, by = 'geneid')


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
library(irrCAC)


rs$common_sig = ifelse(rs$p_common<0.05,1,0)
rs$rare_sig = ifelse(rs$P<0.05,1,0)

top_rank_n = nrow(rs) * prop_causal_genes
rs$common_ranksig = ifelse(rank(rs$p_common)<=top_rank_n,1,0)
rs$rare_ranksig = ifelse(rank(rs$P)<=top_rank_n,1,0)

tb_sig = table(rs$common_sig, rs$rare_sig)
tb_ranksig = table(rs$common_ranksig, rs$rare_ranksig)

#kappa
kappa_sig = as.numeric(kappa_func(as.matrix(tb_sig)))
kappa_rank = as.numeric(kappa_func(as.matrix(tb_ranksig)))

#ac1
ac1_sig = gwet.ac1.table(tb_sig)
ac1_ranksig = gwet.ac1.table(tb_ranksig)

ans = data.frame(c1r1 = tb_sig[2,2], c0r1 = tb_sig[1,2], 
                 c1r0 = tb_sig[2,1], c0r0 = tb_sig[1,1], 
                 kappa_sig = kappa_sig, ac1_sig = ac1_sig$coeff.val,
                 c1r1_rank = tb_ranksig[2,2], c0r1_rank = tb_ranksig[1,2], 
                 c1r0_rank = tb_ranksig[2,1], c0r0_rank = tb_ranksig[1,1], 
                 kappa_rank = kappa_rank, ac1_ranksig = ac1_ranksig$coeff.val,
                 h2 = h2, n = n, prop_causal_genes = prop_causal_genes, 
                 shuffle_rare = shuffle_rare, genetic_model = genetic_model,
                 common_rare_sharing_proportion = common_rare_sharing_proportion,
                 seed_i = seed_i)

write.table(ans, paste0("/data/***/common_rare/simulation/common_magma/ans/h2_",h2,"_geneticmodel_",genetic_model,"_shareprop_",common_rare_sharing_proportion,"_shufflerare_",shuffle_rare,"_n_",n,"_propcausal_",prop_causal_genes,"_seed_",seed_i,".txt"), quote = F, row.names = F, sep='\t')









