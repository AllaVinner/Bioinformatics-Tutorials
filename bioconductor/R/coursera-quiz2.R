
library(AnnotationHub)
BiocManager::install('BSgenome.Hsapiens.UCSC.hg19')
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
library(purrr)

hg <- BSgenome.Hsapiens.UCSC.hg19

###############
# Q1 0.4798807
###############

afreq <-as.list(alphabetFrequency(hg$chr22))
afreq$N <- 0
(afreq$G + afreq$C) / sum(unlist(afreq))


###############
# Q2 0.528866
###############
keys <- c('H3K27me3', 'E003', 'narrowPeak')
ah <- AnnotationHub()
qh <- query(ah, keys)
line <- qh[[1]]

line_chr22 <- keepSeqlevels(line, 'chr22', pruning.mode = 'coarse')
line_chr22
chr22 <- hg$chr22
chr22 <- as(chr22, 'DNAStringSet')
names(chr22) <- 'chr22'
remove(chr22, 'N')
letterFrequency(chr22[line_chr22], 'GC', as.prob = TRUE) %>% mean

###############
# Q3 0.004467924
###############

gcval <- letterFrequency(chr22[line_chr22], 'GC', as.prob = TRUE)
length(gcval)
sigval <- line_chr22$signalValue
cor(gcval, sigval)

###############
# Q4
###############
keys <- c('H3K27me3', 'E003', 'fc.signal')
ah <- AnnotationHub()
qh <- query(ah, keys)
fc <- qh[[1]]
fc_all <- import(fc)
fc_all

fc22 <- import(fc, which=line_chr22)



fc22

fc22$score %>% cor(sigval)

line_chr22
fc_view <- Views(fc22, line_chr22)

fc_view %>% mean
fc_view %>% mean %>% cor(sigval)
ov <- findOverlaps(fc22,line_chr22)
ov
subjectHits(ov)
values(fc22)$region <- subjectHits(ov)

library(dplyr)
values(fc22) %>% as_tibble %>% group_by(region) %>% dplyr::summarize(meanscore = mean(score)) %>% cor(sigval)

values(fc22) %>% as_tibble %>% pull('region') %>% unique

cor(ov$signalValue, ov$score)

?subsetByOverlaps

###############
# Q
###############

sort(line_chr22)

###############
# Q
###############
?letterFrequency



?DNAString
