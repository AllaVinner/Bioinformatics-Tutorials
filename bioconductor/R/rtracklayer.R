library(BiocManager)
BiocManager::install('rtracklayer')
BiocManager::install('Rsamtools')
library(rtracklayer)

?import

# BED files are very similar to granges
# Wig files :: Some sort of signal along the genome
# BigWig is a compressed verstion

ah <- AnnotationHub()

table(ah$rdataclass)

ah.bw <- subset(ah, rdataclass == "BigWigFile" & species == 'Homo sapiens')

bw <- ah.bw[[1]]

# 

# Reads all data
import(bw)

chr22 <- import(bw, which=GRanges("chr22", ranges = IRanges(1, 10^8)))

# lift over
# lets you convert between different versions

# Chain file
# Contain info about converting one specific genome to another specific genome
# "ChainFile"







