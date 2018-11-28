# Abundance trends based on SiteYear population indices
source('01_data_prep.R')

library(lme4)
library(broom)
library(merTools)

allpops <- readRDS("allpops.5.rds")

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
            obsSY10 = length(unique(SiteID[which(obsYearbySite >= 10)]))) %>% 
  filter(posSY5 >= 2 & posYear >= 10)

# write.csv(test, "species_observations.csv", row.names = FALSE)

pops <- allpops %>% 
  filter(CommonName %in% test$CommonName)



pops$YearFact <- as.factor(as.character(pops$Year))
pops$Year <- as.numeric(as.character(pops$Year))
pops$SiteID <- as.factor(as.character(pops$SiteID))
pops$zyear <- pops$Year - 2006


CollInd <- function(temp){
  temp <- temp %>% droplevels()
  mod <- glm(round(Index) ~ YearFact + SiteID - 1, 
             data = temp, family = poisson(link = "log"))
  out <- data.frame(Year = levels(temp$YearFact), 
                    Index = coef(mod)[1:length(levels(temp$YearFact))],
                    zIndex = scale(coef(mod)[1:length(levels(temp$YearFact))])[,1])
}

# CollIndGLMER <- function(temp){
#   temp <- temp %>% droplevels()
#   mod <- glmer(round(Index) ~ zyear +
#                  (1 + zyear|SiteID) +
#                  (1 | YearFact), family = poisson(link = "log"),
#                data = temp)
#   newdat <- temp %>% ungroup() %>% dplyr::select(zyear, YearFact) %>% distinct() %>% mutate(SiteID = "000")
#   out <- list()
#   out$Collated_Index <- data.frame(zyear = newdat$zyear, YearFact = newdat$YearFact, 
#                                    Collated_Index = predict(mod, newdat, re.form = ~ 1 | YearFact))
#   out$mod <- tidy(mod)
#   out$sites <- ranef(mod)$SiteID %>% 
#     mutate(SiteID = row.names(.),
#            CommonName = temp$CommonName[1])
#   
#   tempdf <- data.frame(zyear = seq(min(temp$zyear), max(temp$zyear), length.out = 100),
#                        Year = seq(min(temp$Year), max(temp$Year), length.out = 100),
#                        SiteID = "000", YearFact = "0000")
#   exampPreds <- predictInterval(mod, newdata = tempdf, 
#                                 type = "linear.prediction",
#                                 include.resid.var = FALSE, level = 0.95)
#   tempdf <- cbind(tempdf, exampPreds)
#   out$confint <- tempdf
#   
#   data <- arrange(tempdf, zyear)
#   last <- nrow(data)
#   out$perctrend <- (exp(data$fit[last]) - exp(data$fit[1]))/exp(data$fit[1])
#   return(out)
# }

# TODO: what to do about model fitting with error/convergence issue?
# Hoary Edge Skipper has perfectly correlated random effects, could be removed...?
CollIndGLMER <- function(temp){
  temp <- temp %>% droplevels()
  print(temp$CommonName[1])
  mod <- glmer(round(Index) ~ zyear + meanLL +
                 # zyear * firstyearsurv +
                 (1 |SiteID) +
                 # (1 | YearFact) +
                 (1 | SiteYear), 
               offset = log(TotalModelTime),
               family = poisson(link = "log"),
               data = temp)
  newdat <- temp %>%
    ungroup() %>% 
    dplyr::select(zyear, YearFact) %>% 
    distinct() %>% 
    mutate(SiteID = "000",
           meanLL = mean(temp$meanLL, na.rm = TRUE))
  out <- list()
  out$Collated_Index <- data.frame(zyear = newdat$zyear, 
                                   YearFact = newdat$YearFact, 
                                   meanLL = newdat$meanLL,
                                   Collated_Index = predict(mod, newdat, re.form = ~ 1 | YearFact))
  out$mod <- tidy(mod)
  out$sites <- ranef(mod)$SiteID %>% 
    mutate(SiteID = row.names(.),
           CommonName = temp$CommonName[1])
  
  tempdf <- data.frame(zyear = seq(min(temp$zyear), max(temp$zyear), length.out = 100),
                       Year = seq(min(temp$Year), max(temp$Year), length.out = 100),
                       SiteID = "000", YearFact = "000", meanLL = newdat$meanLL[1], SiteYear = "000")
  exampPreds <- predictInterval(mod, newdata = tempdf, 
                                type = "linear.prediction",
                                include.resid.var = FALSE, level = 0.95)
  tempdf <- cbind(tempdf, exampPreds)
  out$confint <- tempdf
  
  data <- arrange(tempdf, zyear)
  last <- nrow(data)
  out$perctrend <- ((exp(data$fit[last]) / exp(data$fit[1])) ^ (1/(data$Year[last]-data$Year[1])) - 1) * 100
  return(out)
}

effort <- pops %>% 
  ungroup() %>% 
  filter(method == "ukbms") %>% 
  dplyr::select(SiteID, Year, meanLL, meanTime, SurvPerYear) %>% 
  distinct()


popmod <- pops %>%  
  filter(CommonName == "Hoary Edge Skipper") %>%
  # filter(Year != 1995) %>%
  filter(method == "gampred") %>% 
  group_by(CommonName, SiteID) %>% 
  mutate(firstyearsurv = min(unique(zyear)),
         yrpersite = length(unique(Year[which(YearTotal > 0)])),
         yrsincestart = Year - min(Year)) %>% 
  filter(yrpersite >= 5) %>% # this choice makes a difference, including marginal sites lowers the trend
  group_by(CommonName, Year) %>% 
  mutate(siteperyr = length(unique(SiteID[which(YearTotal > 0)]))) %>% 
  # filter(siteperyr >= 3) %>% # increasing this would throw out localized populations
  group_by(CommonName, SiteID, Year) %>% 
  mutate(TotalModelTime = 30 * meanTime,
         TotalSurvTime = SurvPerYear * meanTime) %>% 
  group_by(CommonName) %>% 
  mutate(uniqyr = length(unique(Year)),
         uniqsite = length(unique(SiteID))) %>%
  # filter(uniqyr >= 5,
         # uniqsite >= 3) %>%
  do(results = CollIndGLMER(.)) 

# # kitchen sink
# # gives reasonable species trends that don't vary by site
# # this model should be run in framework that allows 
# # phylogeny and traits to explain differences
# temp <- popmod
# temp$rowid <- paste(temp$CommonName, temp$SiteYear, sep = "_")
# mod <- glmer(round(Index) ~ zyear + meanLL +
#                # zyear * firstyearsurv +
#                (1 + zyear + meanLL|CommonName) +
#                (1 | SiteID:CommonName) +
#                (1 | YearFact:CommonName) +
#                # (1 | YearFact) +
#                # (1 | SiteYear) +
#                (1 | rowid), 
#              offset = log(TotalModelTime),
#              family = poisson(link = "log"),
#              data = temp)


# if not run through CollInd functions, can plot popindex values
plt <- ggplot(popmod, aes(x = Year, y = log(Index), group = SiteID)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~SiteID, scales = "free_y")
plt


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
  
  pc <- data.frame(perctrend = popmod[i, ]$results[[1]]$perctrend)
  pc$CommonName <- popmod[i,]$CommonName
  perclist[[i]] <- pc
}
cidf <- bind_rows(cilist) %>% 
  mutate(Year = as.numeric(as.character(YearFact)))
moddf <- bind_rows(modlist)
sitesdf <- bind_rows(siteslist)
confdf <- bind_rows(confintlist)
percdf <- bind_rows(perclist)

trends <- moddf %>% 
  filter(term == "zyear") %>% 
  arrange(estimate) %>% 
  data.frame()



a <- ggplot(cidf, aes(x = Year, y = exp(Collated_Index)*1800)) +
  geom_point(size = .5) +
  geom_line(data = confdf, aes(x = Year, y = exp(fit) * 30*60), inherit.aes = FALSE, size = .5, linetype = "dashed") +
  geom_ribbon(data = confdf, aes(x = Year, ymin = exp(lwr)*30*60, ymax = exp(upr)*30*60), alpha = .1, inherit.aes = FALSE) +
  expand_limits(y = 0) +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle("Trend of total number of butterflies counted at average site") +
  labs(y = "Predicted number observed") +
  theme_bw(base_size = 8) +
  theme(
    #panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~CommonName, scales = "free_y")

a


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

# abundance model
sitesdf <- ranef(mod)$SiteID %>% 
  mutate(SiteID = row.names(.),
         meanLL = NULL)


sitedat <- sitesdf %>% 
  left_join(sitecov) %>% 
  left_join(sitecov2) %>% 
  # filter(CommonName %in% residents$CommonName) %>% 
  left_join(sitevars) %>% 
  droplevels() %>% 
  ungroup()
names(sitedat)[1] <- "Trend_intc"

# other variables like survey effort not important!
# SiteID random intercept doesn't work with lmer, ==  # levels to nrow
sitemod <- lmer(zyear ~ StartYear + Length + meanLL + meanDur + meanSurv + (1 | SiteID),
                data = sitedat)
summary(sitemod)

sitemod <- lmer(zyear ~ StartYear + Length + (1 | SiteID),
                data = sitedat)
summary(sitemod)

# library(randomForest)
# datrf <- sitedat %>% filter(complete.cases(.))
# sitemod <- randomForest(x = datrf[, -c(2:3)], y = datrf[, 2])


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
  filter(method == "gampred") %>% 
  filter(Year >= 1996,
         is.na(meanTime) == FALSE) %>%
  mutate(zlistlength = meanLL - mean(meanLL, na.rm = TRUE)) %>% 
  # # this uses raw counts with no imputing missing values
  # mutate(TotalTime = meanTime * SurvPerYear) %>% 
  # group_by(SiteID, YearFact, Year, zyear, SiteYear) %>%
  # summarise(Index = sum(YearTotal),
  #           sumTime = mean(TotalTime, na.rm = TRUE),
  #           meanLL = mean(meanLL, na.rm = TRUE))
  # this uses ukbms/gampred indices and imputes missing surveys
  mutate(TotalTime = meanTime * 30) %>% 
    group_by(SiteID, YearFact, Year, zyear, SiteYear) %>%
    summarise(Index = sum(Index),
              sumTime = mean(TotalTime, na.rm = TRUE),
              meanLL = mean(zlistlength, na.rm = TRUE))

temp <- popall %>% droplevels()

# 1st model attempt
mod <- glmer(round(Index) ~ zyear + 
                (1 + zyear|SiteID) +
               (1 | YearFact), 
             offset = log(sumTime),  family = poisson(link = "log"),
             data = temp)
# 2nd attempt with effort covariates
mod <- glmer(round(Index) ~ zyear + meanLL +
               (1 + zyear + meanLL|SiteID) +
               (1 | YearFact) +
               (1 | SiteYear),
             offset = log(sumTime),  family = poisson(link = "log"),
             data = temp)


newdat <- temp %>% 
  ungroup() %>% 
  dplyr::select(zyear, YearFact) %>% 
  distinct() %>% 
  mutate(SiteID = "000",
         meanLL = mean(temp$meanLL, na.rm = TRUE))
Collated_Index <- data.frame(zyear = newdat$zyear, YearFact = newdat$YearFact, 
                                 Collated_Index = predict(mod, newdat, re.form = ~ 1 | YearFact))
Collated_Index$PredTotal <- exp(Collated_Index$Collated_Index) * 30 * 60
Collated_Index$Year <- as.numeric(as.character(Collated_Index$YearFact))

tempdf <- data.frame(zyear = seq(min(temp$zyear), max(temp$zyear), length.out = 100),
                     Year = seq(min(temp$Year), max(temp$Year), length.out = 100),
                     SiteID = "000", YearFact = "0000", meanLL = 0, SiteYear = "0000")
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
  theme(
    #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5))
  
a


poptrendperc <- popmod %>% 
  do(trend = perctrend(.)) %>% 
  unnest()
