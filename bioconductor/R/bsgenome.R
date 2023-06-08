library(BiocManager)
BiocManager::install('BSgenome')
BiocManager::install('BSgenome.Scerevisiae.UCSC.sacCer2')
library(BSgenome)

# Used representing full genomes in bioconductor

available.genomes() # List all available genomes

library(BSgenome.Scerevisiae.UCSC.sacCer2)
Scerevisiae

seqnames(Scerevisiae)
seqlengths(Scerevisiae)
# So far nothing has been loaded into memory

chr1 <- Scerevisiae$chrI

letterFrequency(chr1, 'GC', as.prob = TRUE)

# bsapply
# A bit different than apply (it loades and unloades the genoems as we need them)
param = new('BSParams', X = Scerevisiae, FUN = letterFrequency)
unlist(bsapply(param, 'GC', as.prob = TRUE))




#########
# Views
########

dnaseq <-  DNAString('ACGTACGT')


vi <- matchPattern(dnaseq, Scerevisiae$chrI)
# View
vi

# Materizalize 
ranges(vi) # to ranges

Scerevisiae$chrI[57932: 57939]

# Good to represent a subpart of a bigger object

# functions
shift(vi, 4)

gr <- vmatchPattern(dnaseq, Scerevisiae)

vi2 <- Views(Scerevisiae, gr)
vi2

#########
library(AnnotationHub)

ah <- AnnotationHub()
qh <- query(ah, c('sacCer2', 'genes'))
qh

genes <- ah[["AH7048"]]

proms <- promoters(genes)

# Some genes start right at the beginign, hence the promotors functions
# will give you a negative number
proms <- trim(proms)
proms


promView <- Views(Scerevisiae, proms)

promView

gcProm <- letterFrequency(promView, 'GC', as.prob = TRUE)
gcProm


plot(density(gcProm))



