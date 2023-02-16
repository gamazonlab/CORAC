# -*- coding: utf-8 -*-

import os
import hail as hl
import pandas as pd
import sys

raw_phe = int(sys.argv[1])
ini_num = int(sys.argv[2])
phe_coding = ini_num*100 + raw_phe

mt=hl.read_matrix_table('/data/xxxxxxx_lab/zhoud2/common_rare/results.mt/')

# hail matrixtable to spark Dataframe
data=mt.entries()
data=data.select(data.Pvalue_Burden,data.Pvalue_SKAT,data.BETA_Burden,data.SE_Burden)
data_save=data.to_spark()

#
column_key = pd.read_csv('/data/xxxxxxx_lab/zhoud2/common_rare/phenocode_info.csv', keep_default_na=False)
column_key['filename']=column_key['phenocode']+"_"+column_key['coding']

tmp_data=data_save.filter((data_save["trait_type"]==column_key.loc[2,"trait_type"])
   & (data_save["phenocode"]==column_key.loc[phe_coding,"phenocode"])
   & (data_save["pheno_sex"]==column_key.loc[phe_coding,"pheno_sex"])
   & (data_save["coding"]==column_key.loc[phe_coding,"coding"])
   & (data_save["modifier"]==column_key.loc[phe_coding,"modifier"]))
tmp_data.repartition(1).write.option("header","true").option("sep",",").mode("overwrite").csv("/data/xxxxxxx_lab/zhoud2/common_rare/ukbb_500k/tmp1/"+column_key.loc[phe_coding,"filename"])

