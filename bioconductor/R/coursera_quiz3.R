# Bioconductor


##############
# Q1 5.629627
#############3
library(ALL)
data(ALL)

mean(exprs(ALL[,5]))

##############
# Q2 1045
##############
library(biomaRt)
library("hgu95av2.db")
mart <- useMart(host='https://feb2014.archive.ensembl.org', biomart = "ENSEMBL_MART_ENSEMBL")
ensembl<-useDataset("hsapiens_gene_ensembl",mart)
feature_name<-featureNames(ALL)
listDatasets(mart)
annotation_ALL<-getBM(attributes=c("ensembl_gene_id","affy_hg_u95av2"),filters="affy_hg_u95av2",values=feature_name,mart=ensembl)

summary(annotation_ALL)

ann_summary <- table(annotation_ALL[,2])
sum(ann_summary> 1)

##############
# Q3
##############









