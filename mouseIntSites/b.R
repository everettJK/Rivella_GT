library(gt23)
library(dplyr)
library(RMySQL)
library(GenomicRanges)

dbConn <- dbConnect(MySQL(), group='specimen_management')
sampleData <- dbGetQuery(dbConn, 'select * from gtsp')
dbDisconnect(dbConn)
sampleData <- subset(sampleData, Trial == 'Rivella_NSG_mouse')

intSites <- getDBgenomicFragments(sampleData$SpecimenAccNum, 'specimen_management', 'intsites_miseq') 
intSites <- unlist(GRangesList(lapply(split(intSites, intSites$patient), 
              function(x){ stdIntSiteFragments(x) %>%
                           gt23::collapseReplicatesCalcAbunds() %>%
                           annotateIntSites()
              })))

save(intSites, file = 'data.collapsed.RData')
             
             
             