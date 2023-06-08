library(BiocManager)
BiocManager::install('GenomicFeatures')
BiocManager::install('TxDb.Hsapiens.UCSC.hg19.knownGene')
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

# Genes can have multiple transcripts
# Contains information about genes, transcripts, exons, and coding sequences
# Exons can be part of many transcripts
# Transcripts can have many coding sequences

# Issue when talking about transcritps
# 1. Pre-RNA: Before introns have been cut of
# 2. After splicing

gr <- GRanges(seqnames = "chr1", strand = '+', ranges = IRanges(start = 11874, end = 14409))
gr

# granges with the start and end of the gene (ambigous when a gene can have multiple transcripts)
genes(txdb)
# Exists id (as number) for exon, transcripts, and genes


# Find one gene in the gr
subsetByOverlaps(genes(txdb), gr)
 
# find 3 different transcripts (different exons)
subsetByOverlaps(transcripts(txdb), gr)
subsetByOverlaps(exons(txdb), gr)

# How do I find out how exons create the transcripts

# Here we get the exons for each transcripts
subsetByOverlaps(exonsBy(txdb, by = 'tx'), gr)

# Coding sequence usually mean it is translated to a protein

# Transcripts may have multiple reading frames


subsetByOverlaps(cds(txdb), gr)
# Here we get the coding sequences 

subsetByOverlaps(cdsBy(txdb, by = 'tx'), gr)
# here we get the actual coding seq for the three transcripts
# Only one of the three transcripts have a coding sequence.






