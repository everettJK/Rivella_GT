library(ShortRead)
library(dplyr)

sampleData <- read.csv('sampleData.csv', header = TRUE)
nrow(sampleData) == n_distinct(sampleData$bcSeq)

I1 <- readFastq('/data/sequencing/Illumina-archive/191011_M03249_0020_000000000-CN99T/Data/Intensities/BaseCalls/Undetermined_S0_L001_I1_001.fastq.gz')
R2 <- readFastq('/data/sequencing/Illumina-archive/191011_M03249_0020_000000000-CN99T/Data/Intensities/BaseCalls/Undetermined_S0_L001_R2_001.fastq.gz')

bc.f <- data.frame(sort(table(as.character(I1@sread)), decreasing = TRUE)[1:75])

# Expected barcodes in the top 75 read codes.
sampleData$bcSeq %in% bc.f$Var1 

R2 <- trimTailw(R2, 2, '?', 5)
R2 <- R2[width(R2) >= 100]

# GAAAATC most abundant 1st 7 -- matches expected GAAAATC
R2.1st7 <- narrow(R2, 1, 7)
R2.1st7.tbl <- data.frame(sort(table(as.character(R2.1st7@sread)), decreasing = TRUE)[1:50])

# TCTAGCA most abundant 2nd 7 -- matches match expected TCTAGCA
R2.2nd7 <- narrow(R2, 8, 14)
R2.2nd7.tbl <- data.frame(sort(table(as.character(R2.2nd7@sread)), decreasing = TRUE)[1:50])

R2.1st14 <- narrow(R2, 1, 14)

LTRpassingReads <-  R2[as.character(narrow(R2, 1, 14)@sread) == 'GAAAATCTCTAGCA']
LTRpassingReadsGenomic <- unique(narrow(LTRpassingReads, 15, width(LTRpassingReads))@sread)
names(LTRpassingReadsGenomic) <- paste0('S', 1:length(LTRpassingReadsGenomic))
writeFasta(LTRpassingReadsGenomic, file = 'LTRpassingReadsGenomic.fasta')

