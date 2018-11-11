# Abundance trends based on SiteYear population indices
source('01_data_prep.R')



allpops <- readRDS("allpops.rds")

# 87 species when filtering by >10 sites, >10 years
test <- allpops %>% 
  group_by(CommonName) %>% 
  summarise(meanTot = round(mean(YearTotal)),
            meanIndex = round(mean(Index)),
            posSiteYear = length(which(YearTotal > 0)),
            posSite = length(unique(SiteID[which(YearTotal > 0)])),
            posYear = length(unique(Year[which(YearTotal > 0)]))) %>% 
  filter(posSite > 10, posYear > 10)


pops <- allpops %>% 
  filter(CommonName %in% test$CommonName)



pops$YearFact <- as.factor(as.character(pops$Year))
pops$Year <- as.numeric(as.character(pops$Year))
pops$SiteID <- as.factor(as.character(pops$SiteID))
pops$zyear <- scale(pops$Year)[,1]


CollInd <- function(temp){
  temp <- temp %>% droplevels()
  mod <- glm(round(Index) ~ YearFact + SiteID - 1, 
             data = temp, family = poisson(link = "log"))
  out <- data.frame(Year = levels(temp$YearFact), 
                    Index = coef(mod)[1:length(levels(temp$YearFact))],
                    zIndex = scale(coef(mod)[1:length(levels(temp$YearFact))])[,1])
}


popmod <- pops %>%  
  filter(Year != 1995) %>%
  group_by(CommonName, SiteID) %>% 
  mutate(yrpersite = length(unique(Year[which(YearTotal > 0)])),
         yrsincestart = Year - min(Year)) %>% 
  # filter(yrpersite > 1) %>% 
  group_by(CommonName, Year) %>% 
  mutate(siteperyr = length(unique(SiteID[which(YearTotal > 0)]))) %>% 
  # filter(siteperyr > 2) %>%
  droplevels() %>% 
  group_by(CommonName) %>% 
  mutate(uniqyr = length(unique(Year))) %>%
  filter(uniqyr >= 5) %>%
  do(., CollInd(.)) %>% 
  mutate(Year = as.numeric(as.character(Year)))

poptrend <- popmod %>% 
  group_by(CommonName) %>% 
  do(fit = lm(Index ~ Year, data = .)) %>% 
  broom::tidy(fit) %>% 
  filter(term == "Year") %>% 
  arrange(estimate)

perctrend <- function(data){
  data <- arrange(data, Year)
  pred <-  predict(lm(Index ~ Year, data = data))
  last <- length(pred)
  out <- (exp(pred[last]) - exp(pred[1]))/exp(pred[1])
}

poptrendperc <- popmod %>% 
  group_by(CommonName) %>% 
  do(trend = perctrend(.)) %>% 
  unnest()


a <- ggplot(popmod, aes(x = Year, y = Index)) +
  geom_point()+
  geom_smooth(method = "lm") +
  facet_wrap(~CommonName, scales = "free_y")
a



library(lme4)
# Collated index has problem with zero years, it's not missing data though!
mod <- glmer(round(Index) ~ zyear +
               (1 + zyear|SiteID) +
               (1 | YearFact), family = poisson(link = "log"),
             data = temp)

CollInd <- function(temp){
  temp <- temp %>% droplevels()
  mod <- glmer(round(Index) ~ zyear +
                 (1 + zyear|SiteID) +
                 (1 | YearFact), 
             data = temp, family = poisson(link = "log"))
  out <- data.frame(Year = levels(temp$YearFact), 
                    Index = coef(mod)[1:length(levels(temp$YearFact))],
                    zIndex = scale(coef(mod)[1:length(levels(temp$YearFact))])[,1])
}


# maybe use glmer results to have a few datasets output:
# Year Trend and YearFact variation in "collated index"
# Site Intercept and trend for later use in 




traits <- read.csv("data/speciesphenology.csv", header = TRUE)

poptrend <- poptrend %>% 
  left_join(traits, by = c("species" = "CommonName")) %>% 
  dplyr::select(species:p.value, CombinedLatin:HostCategory, ResStatus, WinterStage) %>% 
  left_join(poptrendperc)
write.csv(poptrend, file = "signifpoptrends.csv", row.names = FALSE)

summary(lm(estimate~BroodsGAMmax + HostCategory + ResStatus + WinterStage, data = poptrend))

# all individuals?

pops <- readRDS("allpops.rds")
pops$YearFact <- as.factor(as.character(pops$Year))
pops$SiteID <- as.factor(as.character(pops$SiteID))
pops$zyear <- scale(pops$Year)[,1]


popmod <- pops %>%  
  filter(Year != 1995) %>%
  filter(species != "Cabbage White") %>%  # with or without Cabbage White
  # group_by(species) %>% 
  # mutate(uniqyr = length(unique(Year))) %>%
  # filter(uniqyr >= 10) %>%
  group_by(SiteID) %>% 
  mutate(yrpersite = length(unique(Year))) %>% 
  group_by(Year) %>% 
  mutate(siteperyr = length(unique(SiteID))) %>% 
  # filter(yrpersite > 1) %>% 
  # siteperyr > 20) %>%
  group_by(SiteID, YearFact, Year) %>%
  summarise(PopIndex = sum(PopIndex)) %>% 
  ungroup() %>% 
  do(., CollInd(.)) %>% 
  mutate(Year = as.numeric(as.character(Year)))
#Pieris = "yes")
# popmod <- rbind(popmod, popmod2)

a <- ggplot(popmod, aes(x = Year, y = exp(Index))) +
  geom_point()+
  geom_smooth(method = "lm") +
  # scale_color_viridis(discrete = TRUE, begin = .2, end = .8) +
  ggtitle("Trend of total # of butterflies counted") +
  labs(y = "Collated Index (butterfly-days)")
a


poptrendperc <- popmod %>% 
  do(trend = perctrend(.)) %>% 
  unnest()
