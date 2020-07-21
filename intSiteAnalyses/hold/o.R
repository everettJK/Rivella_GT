library(tidyverse)

# Read in previously published WAS trial d0 data and rename the subjects and then determin the percentage
# of intSites near oncogenes.

n <- 1
WASintSites_d0 <- readRDS('WAS_d0_intSites.rds')
#WASintSites_d0 <- BostonWAS_m6_PBMC.intSites
WASintSites_d0 <- unlist(GRangesList(lapply(split(WASintSites_d0, WASintSites_d0$patient), 
                                            function(x){ x$patient <- paste('WAS subject', n); n <<- n+1; x})))

intSites <- data.frame(intSites)

wasSitesNearOnco <- 
  data.frame(WASintSites_d0) %>%
  group_by(patient) %>%
  summarise(percentNearOnco = n_distinct(posid[abs(nearestOncoFeatureDist) <= 50000]) / n_distinct(posid)) %>%
  ungroup() %>%
  arrange(percentNearOnco) %>%
  mutate(source = 'WAS') %>%
  data.frame()


humanSitesNearOnco <- 
  group_by(intSites, patient) %>%
  summarise(percentNearOnco = n_distinct(posid[abs(nearestOncoFeatureDist) <= 50000]) / n_distinct(posid)) %>%
  ungroup() %>%
  arrange(percentNearOnco) %>%
  mutate(source = 'Rivella') %>%
  data.frame()




# Create a bar plot comparing the percentage of intSites near oncogenes for human subjects vs WAS d0 subjects.
WASvsHumanSubjects <- 
  bind_rows(humanSitesNearOnco, wasSitesNearOnco) %>%
  arrange(percentNearOnco) %>%
  mutate(patient = factor(patient, levels=unique(patient))) %>%
  ggplot(aes(patient, percentNearOnco, fill=source)) +
  theme_bw() +
  scale_fill_manual(values=c('gray75', 'blue')) +
  geom_bar(stat='identity') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  guides(fill = FALSE) +
  scale_y_continuous(limits = c(0, 0.33), labels = scales::percent) +
  labs(x='Subject', y='Sites near oncogenes') +
  theme(text = element_text(size=16))

