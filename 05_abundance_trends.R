# Abundance trends based on SiteYear population indices
source('01_data_prep.R')

library(lme4)

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

CollIndGLMER <- function(temp){
  temp <- temp %>% droplevels()
  mod <- glmer(round(Index) ~ zyear +
                 (1 + zyear|SiteID) +
                 (1 | YearFact), family = poisson(link = "log"),
               data = temp)
  out <- temp %>% ungroup() %>% dplyr::select(zyear, YearFact) %>% distinct() %>% mutate(SiteID = "000")
  
  out$Index <- predict(mod, out, re.form = ~ 1 | YearFact)
  out$glmer_intc <- fixef(mod)[1]
  out$glmer_slope <- fixef(mod)[2]
  return(out)
}


# TODO:
# BROOM summarise glmer to get slope/se
# YearRE to plot
# SiteRE to plot site covaraites

popmod <- pops %>%  
  filter(Year != 1995) %>%
  group_by(CommonName, SiteID) %>% 
  mutate(yrpersite = length(unique(Year[which(YearTotal > 0)])),
         yrsincestart = Year - min(Year)) %>% 
  filter(yrpersite >= 3) %>%
  group_by(CommonName, Year) %>% 
  mutate(siteperyr = length(unique(SiteID[which(YearTotal > 0)]))) %>% 
  filter(siteperyr >= 3) %>%
  droplevels() %>% 
  group_by(CommonName) %>% 
  mutate(uniqyr = length(unique(Year)),
         uniqsite = length(unique(SiteID))) %>%
  filter(uniqyr >= 5,
         uniqsite >= 5) %>%
  do(., CollIndGLMER(.)) %>% 
  mutate(Year = as.numeric(as.character(YearFact)))

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

CollatedGLMER <- function(temp){
  temp <- temp %>% droplevels()
  m <- glmer(round(Index) ~ zyear +
               (1 + zyear|SiteID) +
               (1 | YearFact), 
             data = temp, family = poisson(link = "log"))
  
  
  
  newdat <- data.frame(zyear = unique(temp$zyear), Index = coef(m)$YearFact[, 1])
  newdat$pred <- predict(m, newdat, type = "response", re.form = NA)
  plt <- ggplot(newdat, aes(x = zyear, y = Index)) +
    geom_point() +
    geom_smooth(method = "lm")
  plt
  
  
  
  
  tempdf <- data.frame(zyear = seq(min(temp$zyear), max(temp$zyear), length.out = 100),
                       SiteID = "000", YearFact = "0000")
  exampPreds <- predictInterval(m, newdata = tempdf, 
                                type = "linear.prediction",
                                include.resid.var = FALSE, level = 0.95)
  tempdf <- cbind(tempdf, exampPreds)
  
  pointdf <- temp %>% ungroup() %>% dplyr::select(zyear, YearFact) %>% distinct() %>% mutate(SiteID = "000")
  pointdf$pred <- predict(m, pointdf, re.form = NA)
  pointdf$pred2 <- predict(m, pointdf, re.form = ~ 1 | YearFact)
  
  ggplot(data=tempdf, aes(x = zyear)) +
  geom_line(aes(y = fit)) +
  geom_ribbon(data=tempdf, aes(ymin = lwr, ymax = upr), fill = "blue", alpha = 0.2)
  
  
  plt <- ggplot(pointdf, aes(x = zyear, y = pred2)) + 
    geom_point() +
    geom_smooth(method = "lm")
  
  pred <- pointdf$pred2
  last <- length(pred)
  out <- (exp(pred[last]) - exp(pred[1]))/exp(pred[1])
  
  out
  results <- list()
  results[[1]] <- summary(mod)
  results[[2]] <- coef(mod)
  results[[3]] <- 
    return(mod)
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
