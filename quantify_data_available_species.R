# ON hold until I get new data through 2016
# then rerun GAMs and do population index all at once

#script combining:
#counts filtering
#gam predictions for a collated index
#differs from poptrends.R by accounting for zero counts at sites
# use sites where sum(total) >= 5, YearSeen at site > 1
# otherwise Site's population indices very low and bring down collated index

source('01_data_prep.R')
allcounts <- list()
for(sp in 1:nrow(allspecies)){
  
  species <- allspecies$CommonName[sp]
  pars <- allspecies %>% filter(CommonName == species)
  
  counts <- data %>% 
    filter(CommonName == species) %>% 
    mutate(DOY = yday(SiteDate),
           Year = year(SiteDate))
  
  #get unique surveys, including those where species not counted (but were at one time)
  survs <- surveys %>% 
    # filter(year(SiteDate) %in% unique(counts$Year)) %>% 
    filter(SiteID %in% unique(counts$SiteID)) %>% 
    mutate(Year = year(SiteDate))
  
  #Add zeros to surveys when species not counted during a survey
  test <- left_join(survs, counts, by = c("SiteID", "SiteDate", "Week", "SeqID", "Year"))
  counts <- test[, c("SiteID", "SiteDate", "Week", "SeqID", "Total", "Year")]
  counts$Total <- plyr::mapvalues(counts$Total, from = NA, to = 0)
  counts <- left_join(counts, covdata)
  counts$temperature[which(counts$temperature < 50)] <- NA
  counts$duration[which(counts$duration == 0)] <- NA
  
  # scaling covariates
  # previously, duration and temperature not that important
  # will just use list-length scaled by site and month
  # this accounts for survey variation while controlling for site and season
  counts <- counts %>% 
    group_by(SiteID) %>% 
    mutate(zlistlength = as.numeric(scale(listlength)))
  
  
  # trying to add GDD instead of ordinal date
  counts <- counts %>% 
    left_join(gdd) %>% 
    group_by(SiteID, Year) %>% 
    mutate(SurvPerYear = length(unique(SeqID)),
           YearTotal = sum(Total),
           CommonName = species,
           CombinedLatin = allspecies$CombinedLatin[sp])
  
  # what if cutoff for inclusion is really open?
  # could filter out sites with lower effort later
  # also this better accounts for SiteYear zero counts, 
  # which are important for population trends
  dat <- counts %>% filter(YearTotal >= 0, SurvPerYear >= 5)
  allcounts[[(length(allcounts) + 1)]] <- dat
}

alldat <- bind_rows(allcounts)

# 87 species when filtering by >10 sites, >10 years
test <- alldat %>% 
  ungroup() %>% 
  dplyr::select(CommonName, CombinedLatin, SiteID, Year, YearTotal, SurvPerYear) %>% 
  distinct() %>% 
  group_by(CommonName, SiteID) %>% 
  mutate(posYearbySite = length(unique(Year[which(YearTotal > 0)])),
         obsYearbySite = length(unique(Year[which(YearTotal >= 0)]))) %>% 
  group_by(CommonName, CombinedLatin) %>% 
  summarise(median_Total_SiteYear = round(median(YearTotal)),
            presence_SiteYear = length(which(YearTotal > 0)),
            presence_Site = length(unique(SiteID[which(YearTotal > 0)])),
            presence_Year = length(unique(Year[which(YearTotal > 0)])),
            presence_Site_5years = length(unique(SiteID[which(posYearbySite >= 5)])),
            presence_Site_10years = length(unique(SiteID[which(posYearbySite >= 10)])),
            zeros_too_Site_5years = length(unique(SiteID[which(obsYearbySite >= 5)])),
            zeros_too_Site_10years = length(unique(SiteID[which(obsYearbySite >= 10)])))

write.csv(test, "OH_species_observations.csv", row.names = FALSE)
