library(gt23)
library(dplyr)
library(RMySQL)

dbConn <- dbConnect(MySQL(), group='specimen_management')
sampleData <- dbGetQuery(dbConn, 'select * from gtsp')
dbDisconnect(dbConn)
sampleData <- subset(sampleData, Trial == 'Rivella_NSG_mouse')


calcSampleAbunds <- function (f) 
{
  f <- data.frame(f)
  f$start <- ifelse(as.character(f$strand) == "+", f$start,  f$end)
  f$end <- f$start
  f$sampleWidth <- paste0(f$sampleName, "/", f$width)
  
  dplyr::group_by(f, sampleName, posid) %>% 
    dplyr::mutate(reads = sum(reads)) %>% 
    dplyr::mutate(estAbund = dplyr::n_distinct(sampleWidth)) %>% 
    dplyr::slice(1) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(-sampleWidth) %>% 
  dplyr::group_by(sampleName) %>% 
    dplyr::mutate(relAbund = (estAbund/sum(estAbund)) * 100) %>% 
  dplyr::ungroup() %>% 
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
}

intSites <- getDBgenomicFragments(sampleData$SpecimenAccNum, 'specimen_management', 'intsites_miseq') 
intSites <- unlist(GRangesList(lapply(split(intSites, intSites$patient), 
              function(x){ stdIntSiteFragments(x) %>%
                           calcSampleAbunds() %>%
                           annotateIntSites()
              })))

save(intSites, file = 'data.RData')
intSites <- data.frame(intSites)
             
             
             