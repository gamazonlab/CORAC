
#ml GCC/8.2.0  OpenMPI/3.1.4 Intel/2019.1.144  IntelMPI/2018.4.274 R/3.6.0

#figure
library(data.table)
require(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)



df = as.data.frame(fread('/data/xxxxxxx_lab/zhoud2/common_rare/cor/result/top100_maf05_kappa.txt'))
df$kappa = ifelse(is.na(df$kappa_ans),0,df$kappa_ans)
#df = df[which(df$p_adj_rare_top_100<0.05),]

# df = df[df$chi_test!=0,]

#df$or = df$table_c0_r0 * df$table_c1_r1 / df$table_c1_r0 / df$table_c0_r1

#spearman_ans_chisq = cor.test(df$chi_test,df$n_eff,method = "spearman")
spearman_ans_or = cor.test(df$or,df$n_eff,method = "spearman")
spearman_kappa = cor.test(df$kappa,df$n_eff,method = "spearman")


# pdf('/data/xxxxxxx_lab/zhoud2/common_rare/cor/chisq_top100.pdf',width = 6, height = 6)
# ggplot(df, aes(x = n_eff, y = chi_test)) +
#   geom_point(size = 0.8) +
#   ggtitle(paste0('Spearman rho = ',round(spearman_ans_chisq$estimate, 3))) +
#   geom_smooth(method = lm, color = 'red') +
#   ylab('chisq_stat_2x2_table')+
#   xlab('n_effective')
# dev.off()

pdf('/data/xxxxxxx_lab/zhoud2/common_rare/cor/figure/top100_maf05_kappa.pdf',width = 8, height = 6)
ggplot(df, aes(x = n_eff, y = kappa)) +
  geom_point(size = 0.8) +
  ggtitle(paste0('Spearman rho = ',round(spearman_kappa$estimate, 3))) +
  geom_smooth(method = lm, color = 'red') +
  ylab("COmmon variant and RAre variant Convergence signature (CORAC) \ndenoted by Cohen's kappa coefficient")+
  xlab('n_effective')
dev.off()


df[,'label'] = NA

idx_1 = which(df$kappa > 0.12)
idx_2 = which(df$kappa < -0.025 & df$n_eff > 3.8e5)
idx_3 = which(df$phenocode %in% c('50','21001'))

idx = unique(c(idx_1, idx_2, idx_3))
df[idx,'label'] = df[idx,'description']

df[,'color'] = '#000000'
dark_col = brewer.pal(12, "Paired")[c(8,2,4)]  #orange blue green
df[idx_1,'color'] = dark_col[2] 
df[idx_2,'color'] = dark_col[1]
df[idx_3,'color'] = brewer.pal(8, "Dark2")[8]

#detailed figure
pdf('/data/xxxxxxx_lab/zhoud2/common_rare/cor/figure/top100_labeled_maf05_kappa.pdf',width = 8, height = 6)
ggplot(df, aes(x = n_eff, y = kappa_ans, label = label)) +  #
  geom_point(size = 0.8, color = df$color) + 
  #ggtitle(paste0('Spearman rho = ',round(spearman_ans_or$estimate, 3))) +
  geom_smooth(method = lm, col = 'red')+
  
  theme_bw() +  #rm background
  theme(panel.grid =element_blank()) +  #rm grids
  
  geom_text_repel( #non overlapped labels
    size=2.5,
    colour = df$color,
    #nudge_x=5e4, #shift to the right
    segment.alpha = 0.2,  #transparent of segment
    min.segment.length = 1,
    segment.size = 0.2,
    #fontface='italic',
    seed = 1
  )+
  ylab("COmmon variant and RAre variant Convergence signature (CORAC) \ndenoted by Cohen's kappa coefficient")+
  xlab('Effective sample size')+
  xlim(0,5e5)+
  ylim(-0.075,0.21)

dev.off()




#--------------h2------------------
# 
# df = as.data.frame(fread('/data/xxxxxxx_lab/zhoud2/common_rare/cor/result/top100_maf05.txt'))
# colnames(df)[12] = 'chi_test'
# 
# df$or = df$table_c0_r0 * df$table_c1_r1 / df$table_c1_r0 / df$table_c0_r1
# 
# spearman_ans_chisq = cor.test(df$chi_test,df$h2_liability,method = "spearman")
# spearman_ans_or = cor.test(df$or,df$h2_liability,method = "spearman")













