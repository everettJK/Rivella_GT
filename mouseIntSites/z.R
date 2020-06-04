library(gt23)
library(dplyr)
library(RMySQL)
library(GenomicRanges)

dbConn <- dbConnect(MySQL(), group='specimen_management')
sampleData <- dbGetQuery(dbConn, 'select * from gtsp')
dbDisconnect(dbConn)
sampleData <- subset(sampleData, Trial == 'Boston_WAS' & Timepoint == 'm6' & CellType == 'PBMC')


intSites <- getDBgenomicFragments(sampleData$SpecimenAccNum, 'specimen_management', 'intsites_miseq') 
intSites <- unlist(GRangesList(lapply(split(intSites, intSites$patient), 
                                      function(x){ stdIntSiteFragments(x) %>%
                                          gt23::collapseReplicatesCalcAbunds() %>%
                                          annotateIntSites()
                                      })))

save(intSites, file = 'BostonWAS_m6_PBMC.RData')