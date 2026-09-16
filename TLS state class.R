# Script description: 

# !/usr/bin/env R
# -*- encoding: utf-8 -*-
# Time:2023/03/01 15:15:42
# Author : Gu Yu

options("repos"= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror="http://mirrors.nju.edu.cn/bioconductor/")
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)

rm(list = ls())
suppressPackageStartupMessages(library(dplyr));library(data.table);library(magrittr);library(tibble);library(stringr)

setwd('/data/home/yu/Yu_project/project_TLS/1.for_bulk/')

survival_data = read.table('/data/home/yu/Yu_common_data/TCGA/TCGA-LIHC.survival.tsv',sep = '\t',header = T)

library(openxlsx)

TLS_1 = read.csv("/data/home/yu/Yu_project/project_TLS/1.for_bulk/TCGA-LIHC_TLS_split_info.csv",header = T)
TLS_2 = openxlsx::read.xlsx("/data/home/yu/Yu_project/project_TLS/1.for_bulk/TCGA-LIHC_TLS_split_info-2.xlsx")

TLS_1_index = TLS_1 %>% filter(Type == 'YES') %>% select(ID) %>% unlist()
names(TLS_1_index) = NULL
TLS_2_index = TLS_2 %>% filter(!Type == 'None') %>% select(Item) %>% unlist()
names(TLS_2_index) = NULL
common_index = intersect(TLS_1_index,TLS_2_index)


No_TLS_1_index = TLS_1 %>% filter(Type == 'None') %>% select(ID) %>% unlist()
names(No_TLS_1_index) = NULL
No_TLS_2_index = TLS_2 %>% filter(Type == 'None') %>% select(Item) %>% unlist()
names(No_TLS_2_index) = NULL
common_index_no = intersect(No_TLS_1_index,No_TLS_2_index)

TLS_patient = survival_data[is.element(survival_data$sample,common_index),c(1,2,4)]
TLS_patient$group = 'TLS'
no_TLS_patient = survival_data[is.element(survival_data$sample,common_index_no),c(1,2,4)]
no_TLS_patient$group = 'No'
TCGA_TLS_ornot_survival = rbind.data.frame(TLS_patient,no_TLS_patient)

write.csv(rbind.data.frame(TLS_patient,no_TLS_patient) %>% .[,-1],'OS_for_XT.csv', row.names = F)

# ------------------ 临床资料整理
clinical = fread('/data/home/yu/Yu_common_data/TCGA/TCGA-LIHC.GDC_phenotype.tsv.gz',sep = '\t',data.table = F)
clinical[1,] %>% t() -> look

want_parameter = c('submitter_id.samples','pathologic_M','pathologic_N','pathologic_T','tumor_stage.diagnoses','neoplasm_histologic_grade')

want_clinical1 = clinical[is.element(clinical$submitter_id.samples,common_index),want_parameter]
want_clinical1$group = 'TLS'
want_clinical2 = clinical[is.element(clinical$submitter_id.samples,common_index_no),want_parameter]
want_clinical2$group = 'No'
want_clinical = rbind.data.frame(want_clinical1,want_clinical2)

Clinical_TLS_ornot = TCGA_TLS_ornot_survival %>% left_join(.,want_clinical,by = c('sample' = 'submitter_id.samples'))

# ------------------------------
load('/data/home/yu/Yu_common_data/TCGA/TCGA-LIHC_gdc_count_2022.Rdata')

Counts = exp
c_rn = Counts$gene_name
Count_TLS_ornot = Counts[,is.element(colnames(Counts),Clinical_TLS_ornot$sample)]
Count_TLS_ornot$gene_id = c_rn

tpm = fread('/data/home/yu/Yu_common_data/TCGA/TCGA-LIHC-tpm-2022.txt',sep = '\t',data.table = F,nThread = 30)
colnames(tpm) %<>% str_sub(.,1,16)
t_rn = tpm$gene_id
Tpm_TLS_ornot = tpm[,is.element(colnames(tpm),Clinical_TLS_ornot$sample)]
Tpm_TLS_ornot$gene_id = t_rn

save('survival_data','clinical','Clinical_TLS_ornot','Count_TLS_ornot','tpm','Tpm_TLS_ornot','Counts',
        file = 'TLS_ornot_bulk_data.rdata')
