# data prep, source this first

# packages to load
library(mclust)
library(lubridate)
library(mgcv)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(viridis)

theme_set(theme_bw(base_size = 14)) 


# for parallel model fitting
# need different packages for windows computers
library(foreach) # for parallelized loops
library(doParallel)

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
  group_by(CommonName, CombinedLatin) %>% 
  summarise(n = sum(Total)) %>% 
  arrange(n)


surveys <- distinct(data[, c("SeqID", "SiteID", "SiteDate", "Week")])



covdata <- data %>%
  group_by(SeqID) %>%
  summarise(listlength = length(which(unique(CommonName) %in% allspecies$CommonName)),
            temperature = mean(c(StartTemp, EndTemp), na.rm = TRUE),
            duration = duration[1]) %>%
  distinct() %>% 
  left_join(surveys)

# impute missing duration values (due to misaligned time entries, etc.)
# one one still is NA because never reported end time for surveys.
covdata <- covdata %>% 
  mutate(Year = year(SiteDate)) %>% 
  group_by(SiteID, Year) %>% 
  mutate(duration = ifelse(is.na(duration) == TRUE | duration == 0,
                           median(duration, na.rm = TRUE),
                           duration))

# site differences
# range of # of species seen at different sites makes me think that 
# list-length should be scaled by month and site to account for survey detection variation
sitelist <- covdata %>% 
  group_by(SiteID) %>% 
  summarise(meanll = mean(listlength),
            sdll = sd(listlength),
            meandur = mean(duration, na.rm = TRUE),
            sddur = sd(duration, na.rm = TRUE),
            nsurvtot = length(unique(SeqID)))




sites <- read.csv("data/OHsites2018update.txt") %>% 
  mutate(SiteID = formatC(Name, width = 3, format = "d", flag = "0"))
gdd <- readRDS("data/dailyDD2016.rds")


gdd <- left_join(gdd, sites) %>% 
  dplyr::select(SiteID, SiteDate, degday530, chill0, lat, lon, maxT, minT) %>% 
  mutate(Year = year(SiteDate),
         DOY = yday(SiteDate)) %>% 
  group_by(SiteID, Year) %>% 
  arrange(DOY) %>% 
  mutate(AccumDD = cumsum(degday530),
         AccumChill = cumsum(chill0))

siteGDD <- gdd %>%
  group_by(SiteID, lat, lon) %>% 
  # mutate(season = case_when(month(SiteDate) %in% c(12, 1, 2) ~ "winter",
  #                           month(SiteDate) %in% c(3, 4, 5) ~ "spring",
  #                           month(SiteDate) %in% c(6, 7, 8) ~ "summer",
  #                           month(SiteDate) %in% c(9, 10, 11) ~ "fall"),
  #        )
  filter(DOY == 365) %>%
  summarise(meanGDD = mean(AccumDD),
            meanChill = mean(AccumChill))

# many ways to cluster sites, but using lat/lon is simplest 
# wanted 4 regions for plotting simplicity
sitemod <- densityMclust(siteGDD[,c(2:3)], G = 1:4, modelNames = "EVV")
siteGDD$region <- as.character(sitemod$classification)
# visualize regional clusters
a <- ggplot(data = siteGDD, aes(x = lon, y = lat, group = region, color = region)) + geom_point()
a

siteGDD$region <- plyr::mapvalues(siteGDD$region, from = c("1", "2", "3", "4"), 
                                  to = c("NE", "NW", "CN", "SW"))
gdd <- gdd %>% 
  left_join(siteGDD[, c("SiteID", "region")])

