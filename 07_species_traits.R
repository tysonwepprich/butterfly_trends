# species traits that explain trend differences?


# Geographic distributions and MAT/MAP at occurrence records
library(raster)
library(rgbif)
library(dplyr)
library(mapr)

data <- readr::read_csv("data/data.trim.csv") %>% 
  mutate(SiteID = formatC(SiteID.x, width = 3, format = "d", flag = "0"),
         SiteDate = lubridate::ymd(SiteDate))
# Two azure species are not distinguished in the monitoring, but could be separated by phenology later
data$CommonName[which(data$CommonName == "Spring/Summer Azure")] <- "Azures"

data <- data %>% 
  filter(year(SiteDate) >= 1996, year(SiteDate) <= 2016)

# filter unidentified species
allspecies <- data %>% 
  filter(CommonName %in% unique(CommonName)[1:122]) %>% 
  group_by(CommonName, CombinedLatin, Genus, Species) %>% 
  summarise(n = sum(Total)) %>% 
  arrange(n)

source('rangefuncs.R')


splist <- c(allspecies$CombinedLatin[1:5])
keys <- sapply(splist, function(x) name_suggest(x)$key[1], USE.NAMES=FALSE)
occ_search(taxonKey=keys, limit=5)

for (sp in 1:nrow(allspecies)){
  spec <- paste(allspecies$Genus[sp], allspecies$Species[sp], sep = " ")
  sp_key <- name_suggest(spec)
  occ <- occ_search(taxonKey = sp_key[1, "key"], hasCoordinate = TRUE, continent = "north_america",
                    hasGeospatialIssue = FALSE, limit = 1000)
  
  # range map
  x <- map_fetch(taxonKey = sp_key, year = 1988:2018)
  map_ggplot(occ, "usa")
  

  }