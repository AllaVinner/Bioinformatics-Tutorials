library(GenomicRanges)
library(AnnotationHub)
library(BiocManager)
library(rtracklayer)
library(purrr)
library(dplyr)
ah <- AnnotationHub()

###################
# Q1 :: 26641
###################

key <- "CpG Islands"

qhs <- query(ah, key)
qhs$genome
# Let's pick the first one (hg19)
genes <- qhs[[1]]


seqlevels(genes)

autosome_names <- map_chr(1:22, ~paste0('chr', as.character(.x)))
autosomes <- keepSeqlevels(genes, autosome_names, pruning.mode = 'coarse')
seqlevels(autosomes)

length(autosomes) # 26641

####################
# Q2 :: 1031
####################

autosomes %>% keepSeqlevels('chr4', pruning.mode = 'coarse') %>% length 


####################
# Q3 :: 41135164
####################

keys <- c('H3K4me3', 'E003', 'narrowPeak')
ah <- AnnotationHub()
qh <- query(ah, keys)
qh
hist <- qh[[1]]
hist

hist %>%  keepSeqlevels(autosome_names, pruning.mode = 'coarse') %>%
 GenomicRanges::disjoin() %>% width %>% sum

################
# Q4 :: 4.770728
################
keys <- c("E003","H3K27me3", "narrowPeak")
ah <- AnnotationHub()
qh <- query(ah, keys)
hist <- qh[[1]]

value_df <- hist %>% keepSeqlevels(autosome_names, pruning.mode = 'coarse') %>% values
mean(value_df$signalValue)

################
# Q5 :: 10289096
################
keys <- c('H3K4me3', 'E003', 'narrowPeak')
ah <- AnnotationHub()
qh <- query(ah, keys)
g1 <- qh[[1]] %>%  keepSeqlevels(autosome_names, pruning.mode = 'coarse')

keys <- c("E003","H3K27me3", "narrowPeak")
ah <- AnnotationHub()
qh <- query(ah, keys)
g2 <- qh[[1]]  %>%  keepSeqlevels(autosome_names, pruning.mode = 'coarse')

GenomicRanges::intersect(g1, g2) %>% GenomicRanges::disjoin() %>% width %>% sum

#################
# Q6 :: 0.5383644
#################

key <- "CpG Islands"
qhs <- query(ah, key)
islands <- qhs[[1]] %>%  keepSeqlevels(autosome_names, pruning.mode = 'coarse')

bivalent <- GenomicRanges::intersect(g1, g2) 

num_overlapping <- countOverlaps(bivalent, islands) %>% purrr::keep(~.x > 0) %>% length() 
num_overlapping / length(bivalent)

#################
# Q7 :: 0.241688
#################

num_overlapping_bases <- GenomicRanges::intersect(islands, bivalent) %>% GenomicRanges::disjoin() %>% width %>% sum
num_island_bases <- GenomicRanges::disjoin(islands) %>% width %>% sum
num_overlapping_bases / num_island_bases

#################
# Q8 :: 9782086
#################


padded_islands <-islands %>% resize(20000+width(.), fix = 'center')

padded_islands %>% GenomicRanges::intersect(bivalent) %>% width %>% sum

#################
# Q9 :: 0.007047481
#################
genome_length <- seqlengths(islands) %>% sum
island_length <- islands %>% GenomicRanges::disjoin() %>% width %>% sum
island_length / genome_length

#################
# Q10 :: 169.0962
#################
inOut <- matrix(0,2, 2)
inOut[1,1] <- sum(width(GenomicRanges::intersect(islands, bivalent, ignore.strand = TRUE)))
inOut[1,2] <- sum(width(GenomicRanges::setdiff(bivalent,islands, ignore.strand = TRUE)))
inOut[2,1] <- sum(width(GenomicRanges::setdiff(islands, bivalent, ignore.strand = TRUE)))
inOut[2,2] <- genome_length - sum(inOut)
inOut

oddsRatio <- inOut[1,1] * inOut[2,2] /(inOut[1,2] * inOut[2,1])
oddsRatio

