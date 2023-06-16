
#BiocManager::install('GEOquery')
library(GEOquery)

# gene expression
# Repository where the dataset has an id number
# Contain alot of different data
eList <- getGEO('GSE11675')

# a list of different datasets
eList

edata <- eList[[1]]

# Expression set
edata
pData(edata)

# Has raw data and processed data.
# Fastq raw data
# SAM processed data

#Get raw data (you get processed data by default)
elist2 <- getGEOSuppFiles('GSE11675')
elist2 # File with the tarred file






