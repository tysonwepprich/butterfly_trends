# Abundance trends based on SiteYear population indices
source('01_data_prep.R')
source('utils.R')

library(lme4)
library(broom)
library(merTools)
library(rtrim)
# devtools::install_github("jknape/poptrend")
library(poptrend)

allpops <- readRDS("allpops.5.rds")

# 87 species after filtering species by data available
test <- allpops %>% 
  filter(method == "ukbms") %>%
  group_by(CommonName, SiteID) %>% 
  mutate(posYearbySite = length(unique(Year[which(YearTotal > 0)])),
         obsYearbySite = length(unique(Year[which(YearTotal >= 0)]))) %>% 
  filter(posYearbySite >= 3, obsYearbySite >= 5) %>% # this choice makes a difference, including marginal sites lowers mean count
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
pops$SiteYear <- as.factor(pops$SiteYear)
pops$Year <- as.numeric(as.character(pops$Year))
pops$SiteID <- as.factor(as.character(pops$SiteID))
pops$zyear <- pops$Year - 2006
pops$zmeanLL <- pops$meanLL - mean(pops$meanLL)
pops$region <- as.factor(stringr::str_split_fixed(pops$RegYear, pattern = coll("_"), 2 )[, 1])


popmod <- pops %>%
  filter(CommonName == "American Copper") %>%
  filter(method == "gampred") %>% 
  group_by(CommonName, SiteID) %>% 
  mutate(firstyearsurv = min(unique(zyear)),
         pospersite = length(unique(Year[which(YearTotal > 0)])),
         obspersite = length(unique(Year)), # includes zero counts
         yrsincestart = Year - min(Year),
         sitezll = zmeanLL - mean(zmeanLL)) %>% 
  group_by(CommonName, Year) %>% 
  mutate(siteperyr = length(unique(SiteID[which(YearTotal > 0)]))) %>% 
  group_by(CommonName, SiteID, Year) %>% 
  mutate(TotalModelTime = 30 * meanTime,
         TotalSurvTime = SurvPerYear * meanTime,
         ModTotal = Index / meanTime) %>% 
  group_by(CommonName) %>% 
  mutate(uniqyr = length(unique(Year)),
         uniqsite = length(unique(SiteID))) %>%
  do(results = CollIndGLMER(.)) 

# years since establishment effect on species list-length?
plt <- ggplot(popmod, aes(x = yrsincestart, y = sitezll)) +
  geom_point() +
  geom_smooth()
plt

CollIndGLMER <- function(temp){
  temp <- temp %>% 
    droplevels() %>% 
    mutate(ztime = sqrt(meanTime) - mean(sqrt(meanTime)))
  print(temp$CommonName[1])
  mod2 <- glmer(round(Index) ~ zyear + zmeanLL + #log(meanTime) +
                 (1 + zyear | SiteID) + # could have meanLL here to vary by site too
                 (1 | YearFact) +
                 (1 | SiteYear), 
               offset = log(meanTime),
               family = poisson(link = "log"),
               data = temp,
               control=glmerControl(optCtrl=list(maxfun=2e6)))
  
  # for species with low sample, simpler model helps convergence
  if(length(mod@optinfo$conv$lme4) > 0){
    mod <- glmer(round(Index) ~ zyear + zmeanLL +
                   (1 | SiteID) + # could have meanLL here to vary by site too
                   (1 | YearFact) +
                   (1 | SiteYear), 
                 offset = log(TotalModelTime), 
                 family = poisson(link = "log"),
                 data = temp,
                 control=glmerControl(optCtrl=list(maxfun=2e6)))
  }
  
  
  newdat <- temp %>%
    ungroup() %>% 
    dplyr::select(zyear, YearFact) %>% 
    distinct() %>% 
    mutate(SiteID = "000",
           zmeanLL = 0,
           ztime = 0)
  out <- list()
  out$model <- mod
  out$data <- temp
  out$modwarnings <- mod@optinfo$conv$lme4
  out$Collated_Index <- data.frame(zyear = newdat$zyear, 
                                   YearFact = newdat$YearFact, 
                                   zmeanLL = newdat$zmeanLL,
                                   Collated_Index = predict(mod, newdat, re.form = ~ 1 | YearFact))
  out$mod <- tidy(mod)
  out$sites <- ranef(mod)$SiteID %>% 
    mutate(SiteID = row.names(.),
           CommonName = temp$CommonName[1])
  
  tempdf <- data.frame(zyear = seq(min(temp$zyear), max(temp$zyear), length.out = 100),
                       Year = seq(min(temp$Year), max(temp$Year), length.out = 100),
                       SiteID = "000", YearFact = "000", zmeanLL = newdat$zmeanLL[1], ztime = newdat$ztime[1], SiteYear = "000")
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

# # using TRIM 

temp <- data %>% 
  droplevels() %>%
  mutate(
    # timeoffset = log(TotalModelTime),
    count = round(Index)) %>% 
  group_by(SiteID) %>% 
  mutate(site_w = mean(count))
print(temp$CommonName[1])

z1 <- rtrim::trim(count ~ SiteID + Year + TotalModelTime, data=temp, model=3, serialcor=TRUE, overdisp=TRUE, weights = "site_w")
summary(z1)
overall(z1)
plot(overall(z1))
plot(index(z1))

z4 <- trim(count ~ SiteID + Year, data=temp, model=2, changepoints="all",
           stepwise=TRUE, serialcor=FALSE, overdisp=TRUE)
summary(z4)




# trying out poptrends from Knape 2016


data <- popmod %>%
  filter(method == "gampred") %>%
  filter(CommonName == "American Copper") %>% 
  mutate(ztime = sqrt(meanTime) - mean(sqrt(meanTime)),
         Index = round(Index))

trFit = ptrend(Index ~ trend(Year, tempRE = TRUE, type = "smooth") + SiteID,
               data = data, family = poisson(link = "log"))
trFit = ptrend(Index ~ trend(Year, tempRE = TRUE, type = "smooth") + SiteID,
               data = data)

trLin = ptrend(Index ~ trend(Year, tempRE = TRUE, type = "smooth", fx = FALSE, k = 6) +
                 s(SiteID, bs = "re") +
                 # SiteID +
                 zmeanLL +
                 ztime, data = data, family = nb(theta = NULL, link = "log"))

plot(trFit)
plot(trLin)
change(trLin, 1996, 2016)

gmod_nb <- gam(round(Index) ~ s(zyear) + zmeanLL + ztime +
                 # s(zyear, SiteID, bs = "fs", k = 5, m = 1) + # could have meanLL here to vary by site too
                 s(SiteID, bs = "re") +
                 s(YearFact, bs = "re"),
               # s(SiteYear, bs = "re"),
               # offset = log(meanTime), 
               # family = poisson(link = "log"),
               family = nb(theta = NULL, link = "log"),
               method = "REML",
               data = data)

# effort <- pops %>% 
#   ungroup() %>% 
#   filter(method == "ukbms") %>% 
#   dplyr::select(SiteID, Year, meanLL, meanTime, SurvPerYear) %>% 
#   distinct()




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
temp <- popmod %>% filter(CommonName == "Monarch", pospersite >= 10)
plt <- ggplot(temp, aes(x = Year, y = YearTotal, group = SiteID)) +
  # geom_point(color = "blue") +
  # geom_smooth(method = "glm", se = FALSE, method.args = list(family = "poisson"), color = "blue") +
  # geom_point(data = temp, aes(x = Year, y = Index, group = SiteID), color = "red") +
  # geom_smooth(data = temp, aes(x = Year, y = Index, group = SiteID), method = "glm", se = FALSE, method.args = list(family = "poisson"), color = "red") +
  geom_point(data = temp, aes(x = Year, y = ModTotal, group = SiteID), color = "green") +
  geom_smooth(data = temp, aes(x = Year, y = ModTotal, group = SiteID), method = "glm", se = FALSE, method.args = list(family = "poisson"), color = "green") +
  facet_wrap(~SiteID, scales = "free_y")
plt



newdat <- temp %>%
  ungroup() %>% 
  dplyr::select(zyear, YearFact, SiteID, SiteYear) %>% 
  distinct() %>% 
  mutate(
    meanLL = mean(temp$meanLL, na.rm = TRUE))
ci <- data.frame(zyear = newdat$zyear, 
                 YearFact = newdat$YearFact, 
                 SiteID = newdat$SiteID,
                 SiteYear = newdat$SiteYear,
                 meanLL = newdat$meanLL,
                 Collated_Index = predict(mod, newdat, re.form = NULL))
# Collated_Index = predict(mod, newdat, re.form = ~ (1 | YearFact) + (1 + zyear |SiteID)))
plt <- ggplot(ci, aes(x = YearFact, y = exp(Collated_Index)*1800, group = SiteID)) +
  geom_point() +
  geom_smooth(method = "glm", se = FALSE, method.args = list(family = "poisson")) +
  # geom_point(data = temp, aes(x = YearFact, y = PopIndex, group = SiteID), color = "red") +
  facet_wrap(~SiteID, scales = "free_y")
plt

errorlist <- list()
cilist <- list()
modlist <- list()
siteslist <- list()
confintlist <- list()
perclist <- list()
for (i in 1:nrow(popmod)){
  err <- popmod[i,]$results[[1]]$modwarnings
  err$CommonName <- popmod[i,]$CommonName
  errorlist[[i]] <- err
  
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
  # geom_line(data = confdf, aes(x = Year, y = exp(fit) * 30*60), inherit.aes = FALSE, size = .5, linetype = "dashed") +
  # geom_ribbon(data = confdf, aes(x = Year, ymin = exp(lwr)*30*60, ymax = exp(upr)*30*60), alpha = .1, inherit.aes = FALSE) +
  geom_smooth(method = "glm", se = TRUE, method.args = list(family = "poisson")) +
  geom_point(size = .5) +
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


# busy slopegraphs of all species trends
a <- ggplot(cidf, aes(x = Year, y = log(exp(Collated_Index)*1800), group = CommonName)) +
  # geom_point(size = .5) +
  geom_line(data = confdf, aes(x = Year, y = log(exp(fit) * 30*60), group = CommonName), inherit.aes = FALSE, size = .5) +
  # geom_ribbon(data = confdf, aes(x = Year, ymin = exp(lwr)*30*60, ymax = exp(upr)*30*60), alpha = .1, inherit.aes = FALSE) +
  # geom_smooth(method = "glm", se = TRUE, method.args = list(family = "poisson")) +
  expand_limits(y = 0) +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle("Trend of total number of butterflies counted at average site") +
  labs(y = "Predicted number observed") +
  theme_bw(base_size = 8) +
  theme(
    #panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5))
# facet_wrap(~CommonName, scales = "free_y")
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
  group_by(SiteID) %>% # testing ideas about nonrandom site placement at initiation biasing trends downward 
  mutate(firstyearsurv = min(unique(zyear)),
         yrsincestart = Year - min(Year)) %>% 
  ungroup() %>% 
  # # this uses raw counts with no imputing missing values
  # mutate(TotalTime = meanTime * SurvPerYear) %>% 
  # group_by(SiteID, YearFact, Year, zyear, SiteYear) %>%
  # summarise(Index = sum(YearTotal),
  #           sumTime = mean(TotalTime, na.rm = TRUE),
  #           meanLL = mean(meanLL, na.rm = TRUE))
  # this uses ukbms/gampred indices and imputes missing surveys
  mutate(TotalTime = meanTime * 30,
         zyrs = yrsincestart - mean(yrsincestart)) %>% 
  group_by(SiteID, YearFact, Year, zyear, SiteYear, zyrs) %>%
  summarise(Index = round(sum(Index)),
            sumTime = mean(TotalTime, na.rm = TRUE),
            meanLL = mean(zmeanLL, na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(ztime = log(sumTime) - mean(log(sumTime), na.rm = TRUE)) %>% 
  group_by(SiteID) %>% 
  mutate(obspersite = length(unique(Year))) %>% 
  filter(obspersite >= 10)



temp <- popall %>% droplevels()

# 1st model attempt
mod <- glmer(round(Index) ~ zyear + 
               (1 + zyear|SiteID) +
               (1 | YearFact) +
               (1 | SiteYear), 
             offset = log(sumTime),
             family = poisson(link = "log"),
             data = temp)
# 2nd attempt with effort covariates
mod <- glmer(round(Index) ~ zyear + meanLL +
               (1 + zyear + meanLL|SiteID) +
               (1 | YearFact) +
               (1 | SiteYear),
             # offset = log(sumTime),  
             family = poisson(link = "log"),
             data = temp)
# 3rd attempt with site initiation covariates
mod <- glmer(round(Index) ~ zyear + meanLL + zyrs +
               (1 + zyear + meanLL|SiteID) +
               (1 | YearFact) +
               (1 | SiteYear),
             offset = log(sumTime),  family = poisson(link = "log"),
             data = temp,
             control=glmerControl(optCtrl=list(maxfun=2e6)))


# effect of offset on trend???
# including list-length seems to have effect on trend's signficance/precision
# time offset does a little worse than ztime as parameter (with coef around .7)

m1a <- glmer(Index ~ zyear + meanLL +
               (1 + zyear + meanLL|SiteID) +
               (1 | YearFact) +
               (1 | SiteYear), 
             family = poisson(link = "log"),
             data = temp)

m2a <- glmer(Index ~ zyear + meanLL +
               (1 + zyear + meanLL|SiteID) +
               (1 | YearFact) +
               (1 | SiteYear),
             offset = log(sumTime),
             family = poisson(link = "log"),
             data = temp)

# could also had ztime into Site RE
m3a <- glmer(Index ~ zyear + ztime + meanLL +
               (1 + zyear + meanLL + ztime|SiteID) +
               (1 | YearFact) +
               (1 | SiteYear),
             family = poisson(link = "log"),
             data = temp)

AIC(m1a, m2a, m3a)

# OR TRIM
temp <- temp %>% 
  mutate(
    # timeoffset = log(TotalModelTime),
    count = round(Index) / sumTime * 60)
# count = round(Index))
z1 <- rtrim::trim(count ~ SiteID + Year, data=temp, model=3, serialcor=TRUE, overdisp=TRUE)
summary(z1)
plot(overall(z1))
plot(index(z1))

# OR POPTREND
tr1 = ptrend(Index ~ trend(Year, tempRE = TRUE, type = "smooth") + SiteID,
             data = temp)

tr2 = ptrend(Index ~ trend(Year, tempRE = TRUE, type = "smooth") +
               s(SiteID, bs = "re") +
               # SiteID +
               s(meanLL) +
               s(ztime), data = temp)

gm <- gam(Index ~ s(Year) + 
            s(YearFact, bs = "re") +
            SiteID +
            s(meanLL) +
            s(ztime), data = temp, family = poisson)

plot(tr1)
plot(tr2)
change(tr1, 1996, 2016)
change(tr2, 1996, 2016)



eff <- ggplot(temp, aes(x = Year, y = ztime, group = SiteID)) + 
  geom_point() +
  geom_smooth() + 
  facet_wrap(~SiteID, scales = "free_y")
eff



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


merBoot <- bootMer(mod, predict, nsim = 100, re.form = NA)



poptrendperc <- popmod %>% 
  do(trend = perctrend(.)) %>% 
  unnest()
