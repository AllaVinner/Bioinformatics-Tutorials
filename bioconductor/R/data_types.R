



# In bio Conductor There are 3 times of data

# Experimental data (facts table)
# Metadata (sample annotation)
# Annotation (gene annotation)


# Annotations
# Give context to an experiment

# E.g. annotate a gene interval

# There are annotation packages.
# Or you can query online resources.



# Expression Set
# Representing experiments (micro array, ngs)

# Contains an expression matrix (num features x num samples) (genes x samples)
# Phenotype data (Sample annotation)
# Feature data (feature annotation)

# eSet
# Contatins multiple matrices (multiple expression sets (or atleast multiple facts tables.))

library(Biobase)
library(ALL) # Experimental data package
#BiocManager::install('ALL')
#BiocManager::install('hgu95av2.db')
#BiocManager::install('airway')

data(ALL)
ALL # expression set

?ALL

# Let's explore it

exprs(ALL)[1:4, 1:4] # Get expression data
head(sampleNames(ALL))
head(featureNames(ALL))

pData(ALL) %>% head()# Phenotype data
featureData(ALL)
ALL$sex # Gives you access to the pheno type data

# Subsetting gives you back the correct expression set
ALL[1:4, 1:5]

# Let's fille up the feature table
ids <- featureNames(ALL)[1:5]
ids

# Use the hgu95av2.db pkg
library(hgu95av2.db)
as.list(hgu95av2ENTREZID[ids])

# Recommended to use pData, but phenoData exists
phenoData(ALL) # Mostly a historical thing ...

###############################
# Summarized Experiment
library(airway)
data(airway)
airway

# colData instead of pData
colData(airway)
extData(airway)

rownames(airway)


assayNames(airway)

assay(airway, 'counts') # facts table 'counts'

length(rowRanges(airway))

# Ranges gives us the exons
elementsLength(rowRanges(airways))

start(airway)

subsetByOverlaps(airways, grange)






