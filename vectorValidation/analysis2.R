library(ShortRead)
library(tidyverse)
library(optparse)
source('lib.R')

option_list = list(
  make_option(c("--R1"), type="character", default=NULL, help="comma delimited list of R1 fastq files", metavar="character"),
  make_option(c("--R2"), type="character", default=NULL, help="comma delimited list of R2 fastq files", metavar="character"),
  make_option(c("--refGenomeFasta"), type="character", default='data/expected.fasta', help="reference genome FASTA path", metavar="character"),   
  make_option(c("--refGenomeBWA"), type="character", default='data/expected', help="path to ref genome BWA database", metavar="character"),   
  make_option(c("--minVariantPhredScore"), type="integer", default=20, help="minimum PHRED score allowed for called varinats", metavar="character"),
  make_option(c("--bwaPath"), type="character", default='~/ext/bwa', help="path to bwa binary", metavar="character"), 
  make_option(c("--megahitPath"), type="character", default='~/ext/megahit/bin/megahit', help="path to megahit binary", metavar="character"), 
  make_option(c("--minBWAmappingScore"), type="integer", default=30, help="minimum BWA mapping score", metavar="character"), 
  make_option(c("--samtoolsBin"), type="character", default='~/ext/samtools/bin', help="path to samtools bin", metavar="character"), 
  make_option(c("--bcftoolsBin"), type="character", default='~/ext/bcftools/bin',  help="path to bcftools bin", metavar="character"),
  make_option(c("--trimQualCode"), type="character", default='5',  help="Min qual trim code", metavar="character"),
  make_option(c("--minTrimmedReadLength"), type="integer", default=25,  help="Min trimmed read length", metavar="integer"),
  make_option(c("--outputDir"), type="character", default='./output',  help="Min qual trim code", metavar="character"))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

#opt$R1 <- 'data/Rivella_S23_L001_R1_001.fastq.gz'
#opt$R2 <- 'data/Rivella_S23_L001_R2_001.fastq.gz'

opt$R1 <- 'data/GTSP4096-1_S41_R1_001.fastq.gz'
opt$R2 <- 'data/GTSP4096-1_S41_R2_001.fastq.gz'


prepareTrimmedReads <- function(R1, R2, qualCode = '5'){
  R1 <- shortRead2DNAstringSet(trimTailw(R1, 2, qualCode, qualCode))
  R2 <- shortRead2DNAstringSet(trimTailw(R2, 2, qualCode, qualCode))
  R1 <- R1[width(R1) >= opt$minTrimmedReadLength]
  R2 <- R2[width(R2) >= opt$minTrimmedReadLength]
  n  <- intersect(names(R1), names(R2))
  R1 <- R1[names(R1) %in% n]
  R2 <- R2[names(R2) %in% n]
  list(R1, R2)
}


# Quality trim reads and create trimmed FASTA files.
r <- prepareTrimmedReads(readFastq(opt$R1), readFastq(opt$R2), qualCode = opt$trimQualCode)


writeFasta(r[[1]], file = file.path(opt$outputDir, 'R1.trimmed.fasta'))
writeFasta(r[[2]], file = file.path(opt$outputDir, 'R2.trimmed.fasta'))


# Align trimmed reads to the reference genome.
system(paste0(opt$bwaPath, ' mem -M ', opt$refGenomeBWA, ' ', file.path(opt$outputDir, 'R1.trimmed.fasta') , ' ', 
              file.path(opt$outputDir, 'R2.trimmed.fasta'),  ' > ', file.path(opt$outputDir, 'genome.sam')))
system(paste0(opt$samtoolsBin, '/samtools view -S -b ', file.path(opt$outputDir, 'genome.sam'), ' > ', 
              file.path(opt$outputDir, 'genome.bam')))
invisible(file.path(opt$outputDir, 'genome.sam'))

# Remove read pairs with mapping qualities below the provided min score.
system(paste0(opt$samtoolsBin, '/samtools view -q ', opt$minBWAmappingScore, ' -b ', file.path(opt$outputDir, 'genome.bam'), 
              ' > ',  file.path(opt$outputDir, 'genome_filt.bam')))

# Retrieve a list of aligned reads.
alignedReadsIDs <- system(paste0(opt$samtoolsBin, '/samtools view ', file.path(opt$outputDir, 'genome_filt.bam'), ' | cut  -f 1 | uniq'), intern = TRUE)


# Build contigs with reads that aligned to the reference genome in the proper orientation and sufficient mapping scores.
writeFasta(r[[1]][names(r[[1]]) %in% alignedReadsIDs], file = file.path(opt$outputDir, 'R1.trimmed.genomeAligned.fasta'))
writeFasta(r[[2]][names(r[[2]]) %in% alignedReadsIDs], file = file.path(opt$outputDir, 'R2.trimmed.genomeAligned.fasta'))


opt$contigs <- megaHitContigs(file.path(opt$outputDir, 'R1.trimmed.genomeAligned.fasta'), 
                              file.path(opt$outputDir, 'R2.trimmed.genomeAligned.fasta'),
                              workDir = file.path(opt$outputDir, 'megahit'), 
                              megahit.path = opt$megahitPath)

refGenomeLength <- width(readFasta(opt$refGenomeFasta))

# Align contigs to reference genome and only retain those that map well (alignment flag == 0 or 16).
writeFasta(opt$contigs, file = file.path(opt$outputDir, 'contigs.fasta'))
system(paste0(opt$bwaPath, ' mem -M ', opt$refGenomeBWA, ' ',  file.path(opt$outputDir, 'contigs.fasta'), ' > ', file.path(opt$outputDir, 'contigs.sam')))
sam <- readLines(file.path(opt$outputDir, 'contigs.sam'))
samLines <- sam[! grepl('^@', sam)]



if(length(samLines) != 0){
  sam <- subset(read.table(textConnection(samLines), sep = '\t', header = FALSE, fill = TRUE), V2 %in% c(0, 16))
  sam$length <- sapply(as.character(sam$V13), samMD2length)
  opt$contigs <- opt$contigs[names(opt$contigs) %in% sam$V1]
  
  if(length(opt$contigs) > 0){
    d <- group_by(sam, V1) %>% top_n(1, length) %>% dplyr::slice(1) %>% 
      summarise(start = V4, length = length, editDist = as.integer(str_extract(V12, '\\d+'))) %>% ungroup()
    
    opt$contigStartPos   <- d[match(names(opt$contigs), d$V1),]$start
    opt$contigEndPos     <- opt$contigStartPos + d[match(names(opt$contigs), d$V1),]$length
    opt$contigsEditDists <- d[match(names(opt$contigs), d$V1),]$editDist
    names(opt$contigs)   <- paste0(names(opt$contigs), ' [', opt$contigsEditDists, ']')
    
    # Check to make sure that the assembler did not create something much longer than the reference.
    opt$contigs <- opt$contigs[! opt$contigEndPos > refGenomeLength + 50]
    
  } else {
    opt$contigStartPos <- NA
    opt$contigEndPos <- NA
    opt$contigsEditDists <- NA
  }
} else 
{
  opt$contigs <- Biostrings::DNAStringSet()
  opt$contigStartPos <- NA
  opt$contigEndPos <- NA
  opt$contigsEditDists <- NA
}



system(paste0(opt$samtoolsBin, '/samtools sort -o ', file.path(opt$outputDir, 'genome_filt.sorted.bam'), ' ', file.path(opt$outputDir, 'genome_filt.bam')))
system(paste0(opt$samtoolsBin, '/samtools index ', file.path(opt$outputDir, 'genome_filt.sorted.bam')))


# Determine the maximum read depth. Overlapping mates will count a position twice.
system(paste0(opt$samtoolsBin, '/samtools depth -d 0 ', file.path(opt$outputDir, 'genome_filt.sorted.bam'), ' -o ', file.path(opt$outputDir, 'genome_filt.sorted.depth')))
maxReadDepth <- max(read.table(file.path(opt$outputDir, 'genome_filt.sorted.depth'), sep = '\t', header = FALSE)[,3])


# Create pileup data file for determining depth at specific locations.
# (!) mpileup will remove duplicate reads.
system(paste0(opt$samtoolsBin, '/samtools mpileup -A -a -Q 0 -o ', file.path(opt$outputDir, 'genome_filt.sorted.pileup'), ' -d ', maxReadDepth, 
              ' -f ', opt$refGenomeFasta, ' ', file.path(opt$outputDir, 'genome_filt.sorted.bam')))


# Determine the percentage of the reference genome covered in the pileup data.
opt$pileupData <- tryCatch({
  read.table(file.path(opt$outputDir, 'genome_filt.sorted.pileup'), sep = '\t', header = FALSE, quote = '')[,1:5]
}, error = function(e) {
  return(data.frame())
})


if(nrow(opt$pileupData) > 0){
  refGenomeLength <- nchar(as.character(readFasta(opt$refGenomeFasta)@sread))
  opt$refGenomePercentCovered <- nrow(subset(opt$pileupData,  V4 >= 1))  / refGenomeLength
  opt$refGenomePercentCovered_5reads <- nrow(subset(opt$pileupData,  V4 >= 5))  / refGenomeLength
  
  # If pileup data could be created then we can try to call variants.
  # --max-depth 100000 
  
  system(paste0(opt$bcftoolsBin, '/bcftools mpileup -A -Ou -f ', opt$refGenomeFasta, ' -d 10 ',
                file.path(opt$outputDir, 'genome_filt.sorted.bam'), ' |  ', opt$bcftoolsBin,  '/bcftools call -mv -Oz ', 
                ' -o ', file.path(opt$outputDir, 'genome_filt.sorted.vcf.gz')))
  
  
  # Read in the variant table created by bcf tools. 
  # We use tryCatch() here because the table may be empty only containing the header information.
  opt$variantTable <- tryCatch({
    ### system(paste0(opt$bcftoolsBin, "/bcftools filter -i'QUAL>", opt$minVariantPhredScore, "' ", 
    system(paste0(opt$bcftoolsBin, "/bcftools filter -i'QUAL>", opt$minVariantPhredScore, " && DP>10' ", 
                  file.path(opt$outputDir, 'genome_filt.sorted.vcf.gz'), " -O z -o ", file.path(opt$outputDir, 'genome_filt.sorted.filt.vcf.gz')))
    
    system(paste0(opt$bcftoolsBin, '/bcftools index ', file.path(opt$outputDir, 'genome_filt.sorted.filt.vcf.gz')))
    
    x <- read.table(file.path(opt$outputDir, 'genome_filt.sorted.filt.vcf.gz'), sep = '\t', header = FALSE, comment.char = '#')
    names(x) <- c('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'OTHER')
    x <- x[, !is.na(names(x))]
    x
  },  error=function(cond) {
    return(data.frame()) 
  })
  
  if(nrow(opt$variantTable) > 0){
    opt$variantTable <- tryCatch({
      # Here we parse the pileup data to create a more informative alt call for variants.
      x <- bind_rows(lapply(1:nrow(opt$variantTable), function(i){
        x <- opt$variantTable[i,]
        p <- parsePileUpString(subset(opt$pileupData, V2 == x$POS)$V5)
        
        # Expand the variant call to include the different possibilities.
        x <- x[rep(1, length(p)),]
        x$ALT <- names(p)
        x$percentAlt <- p
        x <- subset(x, percentAlt > 0)
        x$reads <- subset(opt$pileupData, V2 == x$POS[1])$V4
        
        dplyr::select(x, -ID, -FILTER, -INFO, -FORMAT, -OTHER)
      }))
      x
    }, error=function(cond) {
      stop('Error parsing variant occurrences.')
    })
    
    
    # Select the major variant for each position.
    opt$variantTableMajor <- dplyr::filter(opt$variantTable,  percentAlt >= 0.5 & ! ALT == 'e') %>%
      dplyr::group_by(POS) %>%
      dplyr::top_n(1, wt = percentAlt) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup() %>%
      dplyr::filter(reads >= 5)
    
  }
} else {
  opt$refGenomePercentCovered <- 0
  opt$refGenomePercentCovered_5reads  <- 0
  
  opt$errorCode <- 5
  opt$errorMessage <- 'No pileup or variant data available.'
}




