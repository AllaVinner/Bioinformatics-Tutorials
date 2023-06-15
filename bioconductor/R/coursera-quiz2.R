library(AnnotationHub)
BiocManager::install('BSgenome.Hsapiens.UCSC.hg19')
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
library(BSgenome)
library(purrr)
library(dplyr)

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
# Q4 0.9149614
###############
keys <- c('H3K27me3', 'E003', 'fc.signal')
ah <- AnnotationHub()
qh <- query(ah, keys)
fc <- qh[[1]]
fc_all <- import(fc, which = line_chr22)
fc22 <- keepSeqlevels(fc_all, 'chr22', pruning.mode = 'coarse')
fc22
line_chr22
ov <- findOverlaps(fc_all, line_chr22)
values(fc_all)$group <- subjectHits(ov)
values(fc_all)$width <- fc_all %>% ranges() %>% width()
values(fc_all) %>% as_tibble() %>%
  group_by(group) %>%
  summarize(meanscore = sum(score*width)/sum(width)) %>%
  pull(meanscore) %>%
  cor(line_chr22$signalValue)

###############
# Q5 10914671
###############
keys <- c('H3K27me3', 'E003', 'fc.signal')
ah <- AnnotationHub()
qh <- query(ah, keys)
fc <- qh[[1]]
fc_all <- import(fc)
fc22 <- keepSeqlevels(fc_all, 'chr22', pruning.mode = 'coarse')
fc22[fc22$score >= 1] %>% ranges() %>% width() %>% sum()


###############
# Q6 1869937
###############
keys <- c('H3K27me3', 'E003', 'fc.signal')
ah <- AnnotationHub()
qh <- query(ah, keys)
e3 <- import(qh[[1]]) %>% keepSeqlevels('chr22', pruning.mode = 'coarse')

keys <- c('H3K27me3', 'E055', 'fc.signal')
ah <- AnnotationHub()
qh <- query(ah, keys)
e5 <- import(qh[[1]]) %>% keepSeqlevels('chr22', pruning.mode = 'coarse')

roi <- GenomicRanges::intersect(e3[e3$score <= 0.5], e5[e5$score >= 2])
roi %>% ranges() %>% width() %>% sum()

###############
# Q7 0.8340929
###############
key <- c("hg19", "CpG Islands")
ah <- AnnotationHub()
qhs <- query(ah, key)
islands <- qhs[[1]]
islands <- keepSeqlevels(islands, 'chr22', pruning.mode = 'coarse')
hg <- BSgenome.Hsapiens.UCSC.hg19
chr22 <- hg$chr22

island_dna <- Views(hg, islands)

cs <- letterFrequency(island_dna, 'C')
gs <- letterFrequency(island_dna, 'G')
cgs <- dinucleotideFrequency(island_dna)%>% as_tibble() %>% pull(CG)
region_len <- island_dna %>% ranges() %>% width()

obsfreq <- cgs
expfreq <- cs/region_len*gs

mean(obsfreq/expfreq)

###############
# Q8 27263
###############
tatabox <- 'TATAAA'
chr22 <- hg$chr22
dnaset <- as(chr22, 'DNAStringSet')
matches <- vmatchPattern(tatabox, dnaset)

dnaset <- as(list(chr22, reverseComplement(chr22)), 'DNAStringSet')
matches <- vmatchPattern(tatabox, dnaset)

matches %>% map_dbl(length) %>% sum()



###############
# Q9 139 WRONG
###############

library(TxDb.Hsapiens.UCSC.hg19.knownGene)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
trans <-keepSeqlevels(transcripts(txdb), 'chr22', pruning.mode = 'coarse')
proms <- promoters(trans, upstream=900, downstream=100)
cdseq <-keepSeqlevels(cds(txdb), 'chr22', pruning.mode = 'coarse')

ct <- subsetByOverlaps(proms, cdseq)
ct


ct_dna <- Views(hg, ct)
ct_dnas <- as(ct_dna, 'DNAStringSet')

matches <- vmatchPattern('TATAAA', ct_dnas)

match_size <- map(matches, length)
match_size %>% keep(~.x >0) %>% length()

###########################################################

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
trans <-keepSeqlevels(transcripts(txdb), 'chr22', pruning.mode = 'coarse')
cdseq <-keepSeqlevels(cds(txdb), 'chr22', pruning.mode = 'coarse')
ct <- subsetByOverlaps(trans, cdseq)
proms <- promoters(ct, upstream=900, downstream=100)

promplus <- proms[1:794]
promminus <- proms[795:length(proms)]

ct <- subsetByOverlaps(proms, cdseq)
ct


dnaa <- Views(hg, proms)
dnab <- as(dnaa, 'DNAStringSet')
matches <- vmatchPattern('TATAAA', dnab)

plus_dna <- Views(hg, promplus)
minus_dna <- Views(hg, promminus)

plus_dnas <- as(proms, 'DNAStringSet')

matches <- vmatchPattern('TATAAA', plus_dnas)

match_size <- map(matches, length)
match_size %>% keep(~.x >0) %>% length()



#############################
# Q10 (30 000 something, did not get it correct)
##############################

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
trans <-keepSeqlevels(transcripts(txdb), 'chr22', pruning.mode = 'coarse')
cdseq <-keepSeqlevels(cds(txdb), 'chr22', pruning.mode = 'coarse')
ct <- subsetByOverlaps(trans, cdseq)
proms <- promoters(ct, upstream=900, downstream=100)

ov <- findOverlaps(proms, proms)

for (i in seq_along(ov)){
  q <- queryHits(ov)[[i]]
  s <- subjectHits(ov)[[i]]
  if (q == s) {

  }
  lll <- subsetByOverlaps(proms[q], proms[s])
  print(lll)
}

tibble(a = queryHits(findOverlaps(proms, proms)) ) %>% group_by(a) %>% count() %>% pull(n) %>% sum()
length(a)



