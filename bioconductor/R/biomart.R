#BiocManager::install('biomaRt')
library(biomaRt)

# Choose a data base (marts)
# Then choose the dataset inside the dataset
head(listMarts())

mart <- useMart('ensembl')
# Ensemble is a database you could brows and slowly build up you data query
head(listDatasets(mart))
ens <- useDataset('hsapiens_gene_ensembl', mart)

values <- c('202763_at', '209310_2_at', '207500_at')

# Attributes are what you want to retrives
# Filters filters down the dataa
getBM(attributes = c('ensembl_gene_id', 'affy_hg_u133_plus_2'),
      filters = c('affy_hg_u133_plus_2'), values = values, mart = ens)

# Borws your options
listAttributes(ens)
listFilters(ens)

# Attributes are organized into pages
# You cannot make queries with attributes from different pages.
attributePages(ens)


