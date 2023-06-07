
# Use case:
# Study histon marks, that is set to mark active promoters,
# If truem, then this hison mark should be found at many genes promoters.
# We will study this by lookig at a particular cellline. 

library(GenomicRanges)
library(AnnotationHub)
library(BiocManager)
library(rtracklayer)
#BiocManager::install('rtracklayer')
ah <- AnnotationHub()

ah = subset(ah, species == 'Homo sapiens')
qhs <- query(ah, c('H3K4me3', 'Gm12878'))

# What data should we use?
qhs
# Google some of the titles 
# We will use the wgEncode ...
# And one sample from broad peak, and one from narrow

# Get genomic ranges
gr1 <- qhs[[2]]
gr2 <- qhs[[4]]

# Lets look at the peaks
summary(width(gr1)) # Brod
summary(width(gr2)) # Narrow Seem to force the peaks to be of 150 bases

# Lets pick the narrow one
peaks <- gr2

#? Are these peaks enritched in the promoters?

# lets get a reference sequence
qhs <- query(ah, 'RefSeq')
qhs
qhs$genome

# Let's go with the first one ...
genes <- qhs[[1]]
genes

# Get number of transcripts per gene
table(table(genes$name))
# Most genes have a single transcript

# Let's get our promotors
pro <- promoters(genes)

#? Is these promotores 
# Find if the promotors which overlaps peaks

ov <- findOverlaps(pro, peaks)

length(unique(queryHits(ov)))
length(unique(subjectHits(ov))) 
# Num promotors with a peack

length(subsetByOverlaps(peaks, pro, ignore.strand = TRUE)) / length(peaks) 
length(subsetByOverlaps(pro, peaks, ignore.strand = TRUE)) / length(pro) 

# In most cell types, about 50 % of genes are expressed.

# How many bases does the peaks and pro cover
sum(width(reduce(peaks, ignore.strand = TRUE)))
sum(width(reduce(pro, ignore.strand = TRUE)))
# How big is the overlap
sum(width(intersect(peaks, pro, ignore.strand = TRUE)))

# Lets do with the table

inOut = matrix(0, ncol=2, nrow = 2)
colnames(inOut) = c('Pro in', 'Pro out')
rownames(inOut) = c('Peak in', 'Peak out')

inOut[1,1] <- sum(width(intersect(peaks, pro, ignore.strand = TRUE)))
inOut[1,2] <- sum(width(setdiff(peaks, pro, ignore.strand = TRUE)))
inOut[2,1] <- sum(width(setdiff(pro, peaks, ignore.strand = TRUE)))
inOut[2,2] <- 3*10^9 - sum(inOut)

inOut

oddsRatio <- inOut[1,1] * inOut[2,2] /(inOut[1,2] * inOut[2,1])
oddsRatio

# Does this depend alot on our size of genome = 3G

inOut[2,2] <- 0
inOut[2,2] <- 1.5*10^9 - sum(inOut)
oddsRatio <- inOut[1,1] * inOut[2,2] /(inOut[1,2] * inOut[2,1])
oddsRatio


