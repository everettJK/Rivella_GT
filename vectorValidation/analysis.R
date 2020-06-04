library(ShortRead)
expectedSeqFile <- 'data/expected.fasta'

bwa.path      <- '~/ext/bwa'
samtools.path <- '/home/everett/ext/samtools'
bcftools.path <- '/home/everett/ext/bcftools'
megahit.path  <- '/home/everett/ext/megahit'

R1 <- readFastq('data/Rivella_S23_L001_R1_001.fastq.gz')
R2 <- readFastq('data/Rivella_S23_L001_R2_001.fastq.gz')


# Trim reads with a 10nt sliding window approach and trim where 2 or more bases are below Q30 within a window.
R1 <- trimTailw(R1, 2, '>', 5)
R2 <- trimTailw(R2, 2, '>', 5)

# Convert ShortRead objects to DNA string sets, remove lane info for ids
# and sync reads by ids because a number of reads may of been removed by quality trimming.
shortRead2DNAstringSet <- function(x){
  r <- x@sread
  names(r) <- sub('\\s+.+$', '', as.character(x@id))
  r
}

R1 <- shortRead2DNAstringSet(R1)
R2 <- shortRead2DNAstringSet(R2)
n  <- intersect(names(R1), names(R2))
R1 <- R1[names(R1) %in% n]
R2 <- R2[names(R2) %in% n]

# Write out trimmed reads
writeFasta(R1, file = 'output/R1.trimmed.fasta')
writeFasta(R2, file = 'output/R2.trimmed.fasta')

# Create BWA database with expected vector sequence.
system(paste0(bwa.path ,' index -p output/expected -a bwtsw ', expectedSeqFile))

# Align trimmed reads to expected sequence.
system(paste0(bwa.path, ' mem -M output/expected output/R1.trimmed.fasta output/R2.trimmed.fasta > output/alignment.sam'))

# Convert alignment to bam format 
system(paste0(samtools.path, '/bin/samtools view -S -b output/alignment.sam > output/alignment.bam')) 
invisible(file.remove('output/alignment.sam'))

# Filter out alignments with a mapping quality less than 30 and not propery paired
system(paste0(samtools.path, '/bin/samtools view -q 30 -f 0x2 -b -h output/alignment.bam > output/alignment.filt.bam'))

# Sort filtered bam file and index file result.
system(paste0(samtools.path, '/bin/samtools sort -o output/alignment.filt.sorted.bam  output/alignment.filt.bam'))
system(paste0(samtools.path, '/bin/samtools index output/alignment.filt.sorted.bam'))

# Create read pileup data.
# system(paste0(samtools.path, '/bin/samtools mpileup --output  output/alignment.filt.sorted.pileup --max-depth 100000 -f ', expectedSeqFile, ' output/alignment.filt.sorted.bam')) 
            

# Call variants.
system(paste0(bcftools.path , '/bin/bcftools mpileup -A --max-depth 100000 -Ou -f ', expectedSeqFile, ' output/alignment.filt.sorted.bam | ', 
              bcftools.path,  '/bin/bcftools call -mv -Oz -o output/alignment.filt.sorted.vcf.gz'))

# Filter variant calls on PHRED score >= 20 and read depth >= 10.
system(paste0(bcftools.path, "/bin/bcftools filter -i'QUAL>=20 && DP>=10' output/alignment.filt.sorted.vcf.gz -O z -o output/alignment.filt.sorted.vcf.filt.gz"))


# Read in variants
v <- read.table('output/alignment.filt.sorted.vcf.filt.gz', sep = '\t', header = FALSE, comment.char = '#')
names(v) <- c('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT')
v <- v[,!is.na(names(v))]

