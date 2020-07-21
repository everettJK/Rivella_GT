library(gt23)
library(dplyr)
library(ggplot2)
library(RMySQL)
library(GenomicRanges)
source('lib.R')


# Create / load data sets
#-~-~-~-~-~-~-~--~-~-~--~-~-~--~-~-~--~-~-~--~-~-~--~-~-~--~-~-~--~-~-~--~-~-~--~-~-~--~-~-~--~-~-~--~-~-~-

if (! file.exists('compSites.rds')){
  # Retrieve legacy WAS data.
  legacyData <- readRDS('~/projects/gtVISA_dashboard/data/legacyData.rds')
  legacyData$posid <- paste0(legacyData$seqnames, legacyData$strand, ifelse(legacyData$strand == '+', legacyData$start, legacyData$end))
  legacyData <- subset(legacyData, trial == 'WAS_Paris_Cavazzana')
  legacyData$timePoint <- sub('5.9M', 'M5.9', legacyData$timePoint)

  availableLegacyData <- group_by(legacyData, patient, timePoint, cellType) %>%
                         summarise(nSites = n_distinct(posid)) %>%
                         ungroup()

  # Subset early and late WAS timepoints.
  earlyWAS <- bind_rows(subset(legacyData, timePoint == 'D0' & cellType == 'CD34' & patient == 'pFR01'),
                        subset(legacyData, timePoint == 'D0' & cellType == 'CD34' & patient == 'pFR02'),
                        subset(legacyData, timePoint == 'D0' & cellType == 'CD34' & patient == 'pFR03'),
                        subset(legacyData, timePoint == 'D0' & cellType == 'CD34' & patient == 'pFR04'))
  earlyWAS$date <- 'early'

  lateWAS <- bind_rows(subset(legacyData, timePoint == 'M7'   & cellType == 'PBMC' & patient == 'pFR01'),
                       subset(legacyData, timePoint == 'M5.9' & cellType == 'PBMC' & patient == 'pFR02'),
                       subset(legacyData, timePoint == 'M5.2' & cellType == 'PBMC' & patient == 'pFR03'),
                       subset(legacyData, timePoint == 'M5.7' & cellType == 'PBMC' & patient == 'pFR04'))
  lateWAS$date <- 'late'

  WASdata <- bind_rows(earlyWAS, lateWAS)

  # Retrieve published CGD samples
  CGD_samples <- unique(unlist(stringr::str_extract_all(readLines('~/data/publishedGeneTherapySamples/PRJNA527850_CGD_2019'), 'GTSP\\d+')))
  intSites <- gt23::getDBgenomicFragments(CGD_samples, 'specimen_management', 'intsites_miseq')
  intSites$trial <- 'CGD'
  intSites$dataSource <- 'Illumina'
  intSites$posid <- paste0(seqnames(intSites), strand(intSites), ifelse(as.character(strand(intSites)) == '+', start(intSites), end(intSites)))

  # Subset early and late CGD timepoints.
  earlyCGD <- subset(intSites, timePoint == 'D0')
  earlyCGD$date <- 'early'
  lateCGD <- subset(intSites, timePoint == 'M6' & cellType == 'PBMC')
  lateCGD$date <- 'late'
  CGDdata <- c(earlyCGD, lateCGD)

  # Convert GRanges to data frames and back to merged GRanges to handle metadata columns being in different orders.
  data <- makeGRangesFromDataFrame(bind_rows(WASdata, data.frame(CGDdata)), keep.extra.columns = TRUE)

  rm(WASdata, CGDdata, earlyWAS, earlyCGD, lateWAS, lateCGD, legacyData, intSites, availableLegacyData)
  
  # Standardize fragment boundaries, estimates abundances and annotate sites.
  compSites <- stdIntSiteFragments(data[width(data) >= 10]) %>%
               collapseReplicatesCalcAbunds() %>%
               annotateIntSites()
  
  compSites <- subset(compSites, patient != 'pFR04')
  compSites$patient <- sub('pFR01', 'WAS1', compSites$patient)
  compSites$patient <- sub('pFR02', 'WAS2', compSites$patient)
  compSites$patient <- sub('pFR03', 'WAS3', compSites$patient)
  compSites$patient <- sub('pB101', 'CGD1', compSites$patient)
  compSites$patient <- sub('pN102', 'CGD2', compSites$patient)
  compSites$patient <- sub('pN104', 'CGD3', compSites$patient)
  compSites$patient <- sub('pU103', 'CGD4', compSites$patient)
  
  dbConn <- dbConnect(MySQL(), group='specimen_management')
  sampleData <- dbGetQuery(dbConn, 'select * from gtsp')
  dbDisconnect(dbConn)
  
  compSites <- left_join(data.frame(compSites), select(sampleData, SpecimenAccNum, ConcMeasured, volDNAleft, VCN), by = c("GTSP" = "SpecimenAccNum"))
  compSites$inputDNA <- compSites$ConcMeasured * compSites$volDNAleft
  compSites$inputGenomes <- compSites$inputDNA / 0.0066
  compSites$trial <- sub('WAS_Paris_Cavazzana', 'WAS', compSites$trial)
  compSites$vector <- compSites$trial
    
  saveRDS(compSites, file = 'compSites.rds')
} else {
  compSites <- readRDS('compSites.rds')
}


if(! file.exists('sites.rds')){
  dbConn <- dbConnect(MySQL(), group='specimen_management')
  sampleData <- dbGetQuery(dbConn, 'select * from gtsp where Trial = "Rivella_NSG_mouse"')
  dbDisconnect(dbConn)

  intSites <- getDBgenomicFragments(sampleData$SpecimenAccNum, 'specimen_management', 'intsites_miseq') 
  intSites <- unlist(GRangesList(lapply(split(intSites, intSites$patient), 
                                        function(x){ 
                                            stdIntSiteFragments(x) %>%
                                            collapseReplicatesCalcAbunds() %>%
                                            annotateIntSites()
                                        })))
  intSites$trial <- 'Rivella'
  saveRDS(intSites, file = 'sites.rds')
} else {
  intSites <- readRDS('sites.rds')
}


if(! file.exists('vectorSites.rds')){
  s <- c('GTSP2724', 'GTSP2725', 'GTSP2726', 'GTSP2727', 'GTSP2728', 'GTSP2729', 'GTSP2730', 'GTSP2731', 'GTSP2732')
  v <- c('ALS17', 'ALS17', 'ALS17', 'ALS20', 'ALS20', 'ALS20', 'BB305', 'BB305', 'BB305')
  n <- c('ALS17-1', 'ALS17-2', 'ALS17-3', 'ALS20-1', 'ALS20-2', 'ALS20-3', 'BB305-1', 'BB305-2', 'BB305-3')
  vectorSites <- getDBgenomicFragments(s, 'specimen_management', 'intsites_miseq') 
  vectorSites <- unlist(GRangesList(lapply(split(vectorSites, vectorSites$patient), 
                                        function(x){ 
                                          stdIntSiteFragments(x) %>%
                                            collapseReplicatesCalcAbunds() %>%
                                            annotateIntSites()
                                        })))
  vectorSites$trial <- 'Rivella_vectorSites'
  
  dbConn <- dbConnect(MySQL(), group='specimen_management')
  sampleData <- dbGetQuery(dbConn, 'select * from gtsp')
  dbDisconnect(dbConn)
  
  d <- tibble(s = s, v = v, n= n)
  
  d <- left_join(d, select(sampleData, SpecimenAccNum, ConcMeasured, volDNAleft, VCN), by = c("s" = "SpecimenAccNum"))
  d$inputDNA <- d$ConcMeasured * d$volDNAleft
  d$inputGenomes <- d$inputDNA / 0.0066
                 
  vectorSites$vector <- d[match(vectorSites$GTSP, d$s),]$v
  vectorSites$sampleName <- d[match(vectorSites$GTSP, d$s),]$n
  vectorSites$inputGenomes <- d[match(vectorSites$GTSP, d$s),]$inputGenomes
  
  saveRDS(vectorSites, file = 'vectorSites.rds')
} else {
  vectorSites <- readRDS('vectorSites.rds')
}



# VectorSites refers to the d0 experiments designed to look for differences between vectors. 

c <- data.frame(subset(compSites, timePoint == 'D0'))
c$sampleName <- c$patient
d <- data.frame(vectorSites)

vectorSitesTable <- 
  group_by(bind_rows(d, c), sampleName) %>% 
  summarise(percentSitesNearOnco = sprintf("%.2f%%", (n_distinct(posid[abs(nearestOncoFeatureDist) <= 50000]) / n_distinct(posid))*100)) %>% 
  ungroup()

if(! file.exists('vectorSitesGenomicHeatMap.rds') || ! file.exists('vectorSitesEpiGenomicHeatMap.rds') ){
  c <- data.frame(subset(compSites, timePoint == 'D0'))
  d <- data.frame(vectorSites)
  d$patient <- d$sampleName
  
  if(dir.exists('vectorSitesHeatMap')) unlink('vectorSitesHeatMap')
  createGenomicHeatMapData(bind_rows(d, c), outputDir = 'vectorSitesHeatMap')
  sampleOrder <- c("ALS17-1","ALS17-2","ALS17-3","ALS20-1","ALS20-2","ALS20-3","BB305-1","BB305-2","BB305-3", 
                   "WAS1","WAS2","WAS3","CGD1","CGD2","CGD3","CGD4")
  vectorSitesGenomicHeatMap <- createGenomicHeatMapPlot('vectorSitesHeatMap', sampleOrder = sampleOrder, featuresToDrop =  c("NFR.1M", "aCFR.1M", "onco.100k"))
  
  if(dir.exists('vectorSitesEpiHeatMap')) unlink('vectorSitesEpiHeatMap')
  createEpiGenomicHeatMapData(bind_rows(d, c), outputDir = 'vectorSitesEpiHeatMap')
  vectorSitesEpiGenomicHeatMap <- createGenomicHeatMapPlot('vectorSitesEpiHeatMap', sampleOrder = sampleOrder)
  
  saveRDS(vectorSitesGenomicHeatMap, file = 'vectorSitesGenomicHeatMap.rds')
  saveRDS(vectorSitesEpiGenomicHeatMap, file = 'vectorSitesEpiGenomicHeatMap.rds')
} else {
  vectorSitesGenomicHeatMap <- readRDS('vectorSitesGenomicHeatMap.rds')
  vectorSitesEpiGenomicHeatMap <- readRDS('vectorSitesEpiGenomicHeatMap.rds')
}

ggsave(vectorSitesGenomicHeatMap$gen_plot, width = 7, height = 8, units = 'in', file = 'vectorSitesGenomicHeatMap.pdf')
ggsave(vectorSitesEpiGenomicHeatMap$gen_plot, width = 8, height = 6, units = 'in', file = 'vectorSitesEpiGenomicHeatMap.pdf')



intSites$date <- 'late'
d <- bind_rows(data.frame(intSites), data.frame(compSites))
d$date <- ifelse(d$date == 'early', 'Day 0', 'Month 6')
d <- subset(d, patient != 'pFR04') # too few sites
d$nearOncoGene <- ifelse(abs(d$nearestOncoFeatureDist) <= 50000, TRUE, FALSE)

# Compare Rivella m6 CD34 to WAS/CGD d0 CD34
m1 <- matrix(c(n_distinct(subset(d, trial != 'Rivella' & nearOncoGene == FALSE & date == 'Day 0')$posid),
               n_distinct(subset(d, trial != 'Rivella' & nearOncoGene == TRUE  & date == 'Day 0')$posid),
               n_distinct(subset(d, trial == 'Rivella' & nearOncoGene == FALSE)$posid),
               n_distinct(subset(d, trial == 'Rivella' & nearOncoGene == TRUE)$posid)),
               byrow = TRUE, ncol = 2, dimnames = list(c('WAS/CGD', 'Rivella'), c('notNearOnco', 'nearOnco')))

fisher.test(m1)$p.value

# WAS/CGD near oncoGenes
sprintf("%.2f%%", (m1[1,2]/sum(m1[1,]))*100)

# Rivella near oncoGenes
sprintf("%.2f%%", (m1[2,2]/sum(m1[2,]))*100)

m2 <- matrix(c(n_distinct(subset(d, trial != 'Rivella' & nearOncoGene == FALSE & date == 'Month 6')$posid),
               n_distinct(subset(d, trial != 'Rivella' & nearOncoGene == TRUE  & date == 'Month 6')$posid),
               n_distinct(subset(d, trial == 'Rivella' & nearOncoGene == FALSE)$posid),
               n_distinct(subset(d, trial == 'Rivella' & nearOncoGene == TRUE)$posid)),
             byrow = TRUE, ncol = 2, dimnames = list(c('WAS/CGD', 'Rivella'), c('notNearOnco', 'nearOnco')))

fisher.test(m2)$p.value


pd <- group_by(d, trial, patient, date) %>%
      summarise(nSites = n_distinct(posid), pNearOncoGenes = n_distinct(posid[nearOncoGene == TRUE]) / n_distinct(posid)) %>%
      ungroup() 

dataVis1 <- 
  ggplot(pd, aes(date, pNearOncoGenes, color = trial, group = patient)) +
  theme_bw() +
  geom_point(size = 2.5) +
  geom_line() +
  scale_color_manual(name = 'Group', values = c('red', 'dodgerblue', 'green4')) +
  labs(x = '', y = 'Percent sites near oncogenes') +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0.1, 0.4)) +
  theme(text = element_text(size=16),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

ggsave(dataVis1, file = 'WAS_CGD_Rivella_percentOnco_comp.pdf')


sampleOrder <- c("p0018", "p0036","p2091","p2092","p2094","p2095","p2096","p2882","p2883","p2885","WAS1","WAS2","WAS3","CGD1","CGD2","CGD3","CGD4")

if(! file.exists('lateSitesGenomicHeatMap.rds')){
  createGenomicHeatMapData(subset(d, date == 'Month 6'), outputDir = 'lateSitesHeatMap')
  lateSitesGenomicHeatMap <- createGenomicHeatMapPlot('lateSitesHeatMap', sampleOrder = sampleOrder, featuresToDrop =  c("NFR.1M", "aCFR.1M", "onco.100k"))
  saveRDS(lateSitesGenomicHeatMap, file = 'lateSitesGenomicHeatMap.rds')
} else {
  lateSitesGenomicHeatMap <- readRDS('lateSitesGenomicHeatMap.rds')
}

if(! file.exists('lateSitesEpiGenomicHeatMap.rds')){
  createEpiGenomicHeatMapData(subset(d, date == 'Month 6'), outputDir = 'lateSitesEpiHeatMap')
  lateSitesEpiGenomicHeatMap <- createGenomicHeatMapPlot('lateSitesEpiHeatMap', sampleOrder = sampleOrder)
  saveRDS(lateSitesEpiGenomicHeatMap, file = 'lateSitesEpiGenomicHeatMap.rds')
} else {
  lateSitesEpiGenomicHeatMap <- readRDS('lateSitesEpiGenomicHeatMap.rds')
}


ggsave(lateSitesGenomicHeatMap$gen_plot, width = 7, height = 8, units = 'in', file = 'lateSitesGenomicHeatMap.pdf')
ggsave(lateSitesEpiGenomicHeatMap$gen_plot, width = 8, height = 6, units = 'in', file = 'lateSitesEpiGenomicHeatMap.pdf')




d <- bind_rows(data.frame(intSites), data.frame(compSites))
d <- subset(d, date == 'late')

group_by(d, patient) %>% 
  summarise(percentSitesNearOnco = sprintf("%.2f%%", (n_distinct(posid[abs(nearestOncoFeatureDist) <= 50000]) / n_distinct(posid))*100)) %>% 
  ungroup()
