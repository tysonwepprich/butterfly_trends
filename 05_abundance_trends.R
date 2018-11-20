# Abundance trends based on SiteYear population indices
source('01_data_prep.R')

library(lme4)
library(broom)
library(merTools)

allpops <- readRDS("allpops.rds")

# 87 species when filtering by >10 sites, >10 years
test <- allpops %>% 
  group_by(CommonName, SiteID) %>% 
  mutate(posYearbySite = length(unique(Year[which(YearTotal > 0)])),
         obsYearbySite = length(unique(Year[which(YearTotal >= 0)]))) %>% 
  group_by(CommonName) %>% 
  summarise(medTot = round(median(YearTotal)),
            medIndex = round(median(Index)),
            posSiteYear = length(which(YearTotal > 0)),
            posSite = length(unique(SiteID[which(YearTotal > 0)])),
            posYear = length(unique(Year[which(YearTotal > 0)])),
            posSY5 = length(unique(SiteID[which(posYearbySite >= 5)])),
            posSY10 = length(unique(SiteID[which(posYearbySite >= 10)])),
            obsSY5 = length(unique(SiteID[which(obsYearbySite >= 5)])),
            obsSY10 = length(unique(SiteID[which(obsYearbySite >= 10)]))) #%>% 
  # filter(posSite > 10, posYear > 10)

write.csv(test, "species_observations.csv", row.names = FALSE)

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
  newdat <- temp %>% ungroup() %>% dplyr::select(zyear, YearFact) %>% distinct() %>% mutate(SiteID = "000")
  out <- list()
  out$Collated_Index <- data.frame(zyear = newdat$zyear, YearFact = newdat$YearFact, 
                                   Collated_Index = predict(mod, newdat, re.form = ~ 1 | YearFact))
  out$mod <- tidy(mod)
  out$sites <- ranef(mod)$SiteID %>% 
    mutate(SiteID = row.names(.),
           CommonName = temp$CommonName[1])
  
  tempdf <- data.frame(zyear = seq(min(temp$zyear), max(temp$zyear), length.out = 100),
                       Year = seq(min(temp$Year), max(temp$Year), length.out = 100),
                       SiteID = "000", YearFact = "0000")
  exampPreds <- predictInterval(mod, newdata = tempdf, 
                                type = "linear.prediction",
                                include.resid.var = FALSE, level = 0.95)
  tempdf <- cbind(tempdf, exampPreds)
  out$confint <- tempdf
  
  data <- arrange(tempdf, zyear)
  last <- nrow(data)
  out$perctrend <- (exp(data$fit[last]) - exp(data$fit[1]))/exp(data$fit[1])
  return(out)
}


popmod <- pops %>%  
  filter(Year != 1995) %>%
  group_by(CommonName, SiteID) %>% 
  mutate(yrpersite = length(unique(Year[which(YearTotal > 0)])),
         yrsincestart = Year - min(Year)) %>% 
  filter(yrpersite >= 5) %>%
  group_by(CommonName, Year) %>% 
  mutate(siteperyr = length(unique(SiteID[which(YearTotal > 0)]))) %>% 
  filter(siteperyr >= 5) %>%
  droplevels() %>% 
  group_by(CommonName) %>% 
  mutate(uniqyr = length(unique(Year)),
         uniqsite = length(unique(SiteID))) %>%
  filter(uniqyr >= 5,
         uniqsite >= 5) %>%
  do(results = CollIndGLMER(.)) 

# 
# poptrend <- popmod %>% 
#   group_by(CommonName) %>% 
#   do(fit = lm(Index ~ Year, data = .)) %>% 
#   broom::tidy(fit) %>% 
#   filter(term == "Year") %>% 
#   arrange(estimate)
# 
perctrend <- function(data){
  data <- arrange(data, Year)
  pred <-  predict(lm(Index ~ Year, data = data))
  last <- length(pred)
  out <- (exp(pred[last]) - exp(pred[1]))/exp(pred[1])
}
# 
# poptrendperc <- popmod %>% 
#   group_by(CommonName) %>% 
#   do(trend = perctrend(.)) %>% 
#   unnest()

# 
cilist <- list()
modlist <- list()
siteslist <- list()
confintlist <- list()
perclist <- list()
for (i in 1:nrow(popmod)){
  ci <- popmod[i, ]$results[[1]]$Collated_Index
  ci$CommonName <- popmod[i,]$CommonName
  cilist[[i]] <- ci
  
    ml <- popmod[i, ]$results[[1]]$mod
  ml$CommonName <- popmod[i,]$CommonName
  modlist[[i]] <- ml
  
  si <- popmod[i, ]$results[[1]]$sites
  si$CommonName <- popmod[i,]$CommonName
  siteslist[[i]] <- si
  
  cf <- popmod[i, ]$results[[1]]$confint
  cf$CommonName <- popmod[i,]$CommonName
  confintlist[[i]] <- cf
  
  # pc <- popmod[i, ]$results[[1]]$perctrend
  # pc$CommonName <- popmod[i,]$CommonName
  # perclist[[i]] <- pc
}
cidf <- bind_rows(cilist)
moddf <- bind_rows(modlist)
sitesdf <- bind_rows(siteslist)
confdf <- bind_rows(confintlist)
# percdf <- bind_rows(perclist)

trends <- moddf %>% 
  filter(term == "zyear") %>% 
  arrange(estimate) %>% 
  data.frame()


species <- data %>% 
  dplyr::select(CommonName, Genus, Species) %>% 
  distinct() %>% 
  mutate(Latin = paste(Genus, Species, sep = " ")) 
trends <- left_join(trends, species)
outdf <- trends %>% 
  dplyr::select(CommonName, Latin, estimate, std.error, p.value)

write.csv(outdf, "Trends.csv", row.names = FALSE)



a <- ggplot(popmod, aes(x = Year, y = Index)) +
  geom_point()+
  geom_smooth(method = "lm") +
  facet_wrap(~CommonName, scales = "free_y")
a



  pointdf <- temp %>% ungroup() %>% dplyr::select(zyear, YearFact) %>% distinct() %>% mutate(SiteID = "000")
  pointdf$pred <- predict(m, pointdf, re.form = NA)
  pointdf$pred2 <- predict(m, pointdf, re.form = ~ 1 | YearFact)
  
  ggplot(data=tempdf, aes(x = zyear)) +
    geom_line(aes(y = fit)) +
    geom_ribbon(data=tempdf, aes(ymin = lwr, ymax = upr), fill = "blue", alpha = 0.2)



# Site models
sitecov <- read.csv("data/pest_bfly_buff_overtime_MOA.csv") %>% 
  mutate(SiteID = formatC(site, width = 3, format = "d", flag = "0")) %>% 
  filter(buffer == 2000, MOAuse %in% c("her_G", "ins_4")) %>% 
  group_by(SiteID, MOAuse) %>% 
  summarise(meanuse = mean(kg_buff, na.rm = TRUE)) %>% 
  spread(key = MOAuse, value = meanuse, fill = NA)
sitecov2 <- read.csv("data/LC_bflybuff_overtime.csv") %>% 
  mutate(SiteID = formatC(site, width = 3, format = "d", flag = "0")) %>% 
  filter(buffer == 2000, reclass %in% c(2, 4, 6)) %>% 
  group_by(SiteID, reclass) %>% 
  summarise(meanlanduse = mean(km2, na.rm = TRUE)) %>% 
  spread(key = reclass, value = meanlanduse, fill = NA)
names(sitecov2)[2:4] <- c("developed", "forest", "agriculture") 

traits <- read.csv("data/speciesphenology.csv", header = TRUE)
residents <- traits %>% filter(ResStatus == "resident")

sitevars <- pops %>% 
  ungroup() %>% 
  dplyr::select(SiteID, Year, lat, lon, meanLL, meanTime, SurvPerYear, zyear) %>% 
  distinct() %>% 
  group_by(SiteID, lat, lon) %>% 
  summarise(meanLL = mean(meanLL),
            meanDur = mean(meanTime, na.rm = TRUE),
            meanSurv = mean(SurvPerYear),
            StartYear = min(zyear),
            EndYear = max(zyear),
            Length = max(Year) - min(Year))

sitedat <- sitesdf %>% 
  left_join(sitecov) %>% 
  left_join(sitecov2) %>% 
  filter(CommonName %in% residents$CommonName) %>% 
  left_join(sitevars)
names(sitedat)[1] <- "Trend_intc"

# other variables like survey effort not important!
sitemod <- lmer(zyear ~ StartYear + Length + meanLL + meanSurv + (1 | SiteID),
                data = sitedat)
summary(sitemod)

sitemod <- lmer(zyear ~ StartYear + Length + (1 | SiteID),
                data = sitedat)
summary(sitemod)

sitere <- data.frame(Mean_effect = ranef(sitemod)$SiteID[, 1], SiteID = row.names(ranef(sitemod)$SiteID)) %>% 
  left_join(distinct(sitedat[, c("SiteID", "lat", "lon")]))

GGally::ggpairs(sitedat[, c("agriculture", "her_G", "ins_4")])
  pairs(sitedat[, c("developed", "forest", "agriculture", "her_G", "ins_4")])
  
  us <- c(left = -85.15, bottom = 38.35, right = -80, top = 42.1)

  mapplt <- ggplot(sitere, aes(x = lon, y = lat, color = Mean_effect)) +
    geom_point(size = 6, alpha = .7) +
    scale_color_viridis(name = "Site Effect", begin = 0, end = .9, discrete = FALSE) +
    theme_minimal(base_size = 24) +
    theme(legend.position = c(.88, .15),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
  geom_polygon(data = reg.df, aes(x = long, y = lat, group = group), fill = NA, color = "dark gray", inherit.aes = FALSE, size = 1) +
    coord_fixed(1.3, xlim = c(-85.15, -80.2), ylim = c(38.35, 42.1), expand = FALSE, clip = "on") +
    xlab("Longitude") +
    ylab("Latitude")
mapplt

  poptrend <- poptrend %>% 
  left_join(traits, by = c("species" = "CommonName")) %>% 
  dplyr::select(species:p.value, CombinedLatin:HostCategory, ResStatus, WinterStage) %>% 
  left_join(poptrendperc)
write.csv(poptrend, file = "signifpoptrends.csv", row.names = FALSE)



traits <- read.csv("data/speciesphenology.csv", header = TRUE)
trenddat <- trends %>% 
  left_join(traits)
summary(lm(estimate~BroodsGAMmax + HostCategory + ResStatus + WinterStage, data = trenddat))

# all individuals?


popall <- pops %>%  
  ungroup() %>% 
  filter(Year >= 1996,
         is.na(meanTime) == FALSE) %>%
  # filter(CommonName != "Cabbage White") %>%  # with or without Cabbage White
  mutate(TotalTime = meanTime * SurvPerYear) %>% 
  group_by(SiteID, YearFact, Year, zyear) %>%
  summarise(Index = sum(YearTotal),
            sumTime = mean(TotalTime, na.rm = TRUE)) #%>% 

temp <- popall %>% droplevels()

mod <- glmer(round(Index) ~ zyear + 
                (1 + zyear|SiteID) +
               (1 | YearFact), offset = log(sumTime),  family = poisson(link = "log"),
             data = temp)
newdat <- temp %>% ungroup() %>% dplyr::select(zyear, YearFact) %>% distinct() %>% mutate(SiteID = "000")
Collated_Index <- data.frame(zyear = newdat$zyear, YearFact = newdat$YearFact, 
                                 Collated_Index = predict(mod, newdat, re.form = ~ 1 | YearFact))
Collated_Index$PredTotal <- exp(Collated_Index$Collated_Index) * 30 * 60
Collated_Index$Year <- as.numeric(as.character(Collated_Index$YearFact))

tempdf <- data.frame(zyear = seq(min(temp$zyear), max(temp$zyear), length.out = 100),
                     Year = seq(min(temp$Year), max(temp$Year), length.out = 100),
                     SiteID = "000", YearFact = "0000")
intervals <- predictInterval(mod, newdata = tempdf, 
                              type = "linear.prediction",
                              include.resid.var = FALSE, level = 0.95)
tempdf<- cbind(tempdf, intervals)

a <- ggplot(Collated_Index, aes(x = Year, y = PredTotal)) +
  geom_point(size = 5) +
    geom_line(data = tempdf, aes(x = Year, y = exp(fit) * 30*60), inherit.aes = FALSE, size = 2, linetype = "dashed") +
  geom_ribbon(data = tempdf, aes(x = Year, ymin = exp(lwr)*30*60, ymax = exp(upr)*30*60), alpha = .1, inherit.aes = FALSE) +
  expand_limits(y = 0) +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle("Trend of total number of butterflies counted at average site") +
  labs(y = "Predicted number observed") +
  theme_bw(base_size = 24) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5))
  
a


poptrendperc <- popmod %>% 
  do(trend = perctrend(.)) %>% 
  unnest()
