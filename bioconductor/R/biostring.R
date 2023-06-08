



#Biostring Package
# Package for manipulating biological strings (DNA, RNA, aminoacids)
# DNAString, RNAString, AAString, (BIOString, XString (any))
library(BiocManager)
BiocManager::install('Biostrings')
library(Biostrings)


dna1 <-  DNAString('ACGT-G')
dna2 <- DNAStringSet(c('ACG', 'ACGT', 'ACGTT'))

# Looks like a character string but does some optimizations under the hood.
# The dnaset is a collection of DNAString
# Restricted to have characters in the IUPAC_CODE_MAP
IUPAC_CODE_MAP

dna2[[2]] # Get DNAString from DNAStringSet

# Functions
width(dna2)
sort(dna2)
rev(dna1)
reverse(dna2)
reverseComplement(dna2)
translate(dna2)
alphabetFrequency(dna2) # frequence for each DNAString and every symbole
letterFrequency(dna2, letters = 'GC') # G or C
dinucleotideFrequency(dna2) # Count all 2-merse
consensusMatrix(dna2) 


# Matching

# Short read alignment in your terminal
library(BSgenome)
library(BSgenome.Scerevisiae.UCSC.sacCer2)

dna <- DNAString('ACGTACGT')
# match string, set string with string set string

# string sting matching
matchPattern(dna, Scerevisiae$chrI) # return a view (looks like a dna set)
countPattern(dna, Scerevisiae$chrI)

# string to string set
vmatchPattern(dna, Scerevisiae)
# Considers both + and - strand

matchPWM() # position weight matrix
pairwiseAlignment() # Implement a pairwise alignments (Good for small read and genome)
trimLRPatterns() # Trim off sequence adapters.






