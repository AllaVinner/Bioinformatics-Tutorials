
# Many genes can be seen as intervals.
# Granges is a data structure to deal with these intervals

# Two main packaes
# IRanges
# GenomicRanges

############
# IRanges
############

# Iranges - Vector which stores integer intervlas
# S4 object.

#library(BiocManager)
#BiocManager::install('IRanges')

library(IRanges)

# Construct an IRanges
ir <- IRanges(c(1,3,7), width = 4)
ir2 <- IRanges(c(1,2,8), width = 4)
ir

# Properties
start(ir)
end(ir)
width(ir)
length(ir)
ir[1]


# Functions
reduce(ir)
disjoin(ir)
resize(ir, width = 1, fix = 'center')
union(ir, ir2)
ov <- findOverlaps(ir, ir2)
ov
# Adjacency matrix.
# Index of the overlaps
# Read as query overlap subject. 

queryHits(ov)
countOverlaps(ir, ir2)
nearest(ir, ir2)

#####################
# GenomicRanges
#####################
# IRanges with some extra info
# Seqname - Strand  - Ranges
# Rle - Run length encoding

BiocManager::install('GenomicRanges')
library(GenomicRanges)

gr <- GRanges(seqname = c("ch1"), strand = c('+', '-', '+'),
              ranges = IRanges(start = c(1,3,5), width = 3))
gr

# Functions
flank(gr, 5) # Take the ranges before the irange (+ for +, - for - strand)
promoters(gr) # Take 2000 interval downstream, 200 upstream

# Properties
seqinfo(gr)
seqlengths(gr) <- 10
seqlevels(gr)
seqnames(gr)
gaps(gr) # What is not covvered by the granges
sort(gr) # Sort relative to the seq levels too (Usually by chromosome)
genome(gr) <- 'hg19'

# GRanges handles the bookeeping when things are not compatible.
gr2 = gr
genome(gr2) <- 'hg18'
findOverlaps(gr, gr2)
subsetByOverlaps(gr, gr2)

seqlevels(gr, force = TRUE) = 'chr1' # Drops all other chromosomes
dropSeqlevels(gr, 'chr1')
keepSeqlevels(gr, 'chr2')
keepStandardChromosomes(gr)
# Style of writing chromosomes are differnet ... (chr1, Chr1,...)
newStyle <- mapSeqlevels(seqlevels(gr), 'NCBI')
gr <- renameSeqlevels(gr, newStyle)

################################
# DataFrame can store any type which implements 'length;
############################

# Constructor
ir <- IRanges(start = 1:3, width = 2)
df <- DataFrame(ir = ir, score = rnorm(3))
df

# 
gr <- GRanges(seqnames = 'chr1', strange = c('+', '-', '+'),
              ranges = IRanges(start = c(1,3,5), width = 3))

# The values takes a data frame
values(gr) <- DataFrame(score = rnorm(3))
gr
gr$score

makeGRangesFromDataFrame(df, keep.extra.columns = TRUE)




