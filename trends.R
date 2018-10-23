# ON hold until I get new data through 2016
# then rerun GAMs and do population index all at once

#script combining:
#counts filtering
#gam predictions for a collated index
#differs from poptrends.R by accounting for zero counts at sites
# use sites where sum(total) >= 5, YearSeen at site > 1
# otherwise Site's population indices very low and bring down collated index

# packages to load
library(mclust)
library(lubridate)
library(mgcv)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)

# for parallel model fitting
# need different packages for windows computers
library(foreach) # for parallelized loops
if(.Platform$OS.type == "unix"){
  library(doMC)    # parallel backend for foreach, only for linux/mac
}else if(.Platform$OS.type == "windows"){
  library(doSNOW)
}

data <- readr::read_csv("data/data.trim.csv") %>% 
  mutate(SiteID = formatC(SiteID.x, width = 3, format = "d", flag = "0"),
         SiteDate = lubridate::ymd(SiteDate))
data$CommonName[which(data$CommonName == "Spring/Summer Azure")] <- "Azures"

# filter unidentified species
allspecies <- data %>% 
  filter(CommonName %in% unique(CommonName)[1:122]) %>% 
  group_by(CommonName) %>% 
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




sites <- read.csv("data/OHsites_reconciled_update2016.csv") %>% 
  mutate(SiteID = formatC(Name, width = 3, format = "d", flag = "0"))
gdd <- readRDS("data/dailyDD.rds")


gdd <- left_join(gdd, sites) %>% 
  dplyr::select(SiteID, SiteDate, degday1030, chill0, lat, lon, maxT, minT) %>% 
  mutate(Year = year(SiteDate),
         DOY = yday(SiteDate)) %>% 
  group_by(SiteID, Year) %>% 
  arrange(DOY) %>% 
  mutate(AccumDD = cumsum(degday1030),
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





# Fit a generalized additive model (GAM) for each species
# This model will be used to impute missing counts a la UKBMS
# OR be used to compute site x year population indices



ncores <- 30
if(.Platform$OS.type == "unix"){
  registerDoMC(cores = ncores)
}else if(.Platform$OS.type == "windows"){
  cl <- makeCluster(ncores)
  registerDoSNOW(cl)
}

mcoptions <- list(preschedule = FALSE)

# foreach loop
outfiles <- foreach(sp = 1:nrow(allspecies),
                    .combine='c',
                    .packages= c("mgcv", "dplyr", "tidyr", "purrr",
                                 "lubridate", "mclust", "mixsmsn"),
                    .export = c("data", "surveys", "covdata", "allspecies", "gdd"),
                    .inorder = FALSE,
                    .options.multicore = mcoptions) %dopar% {
                      
                      species <- allspecies$CommonName[sp]
                      
                      counts <- data %>% 
                        filter(CommonName == species) %>% 
                        mutate(DOY = yday(SiteDate),
                               Year = year(SiteDate))
                      
                      #get unique surveys, including those where species not counted (but were at one time)
                      survs <- surveys %>% 
                        filter(year(SiteDate) %in% unique(counts$Year)) %>% 
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
                               YearTotal = sum(Total))
                      
                      
                      # if(cutoff == "strict"){
                      #   datGAM <- counts[YearTotal >= 3]
                      #   datGAM <- datGAM[SurvPerYear >= 15]
                      # }
                      # if(cutoff == "loose"){
                      #   datGAM <- counts %>% filter(YearTotal >= 1, SurvPerYear >= 10)
                      # }
                      
                      # what if cutoff for inclusion is really open?
                      dat <- counts %>% filter(YearTotal >= 1, SurvPerYear >= 5)
                      
                      
                      
                      if(nrow(datGAM) == 0){
                        mod <- NA
                      }else{
                        
                        dat$Year <- as.factor(as.character(dat$Year))
                        dat$region <- as.factor(as.character(dat$region))
                        dat$SiteID <- as.factor(as.character(dat$SiteID))
                        dat$SiteYear <- as.factor(paste(dat$SiteID, dat$Year, sep = "_"))
                        dat$zlistlength[which(is.na(dat$zlistlength))] <- 0
                        dat$RegYear <- as.factor(paste(dat$region, dat$Year, sep = "_"))
                        dat <- as.data.frame(dat)
                        
                        dat <- dat[-which(is.na(dat$AccumDD)), ]
                        
                        # #silly filters for univoltine species with outliers
                        # if(species == "Baltimore"){
                        #   dat <- dat[-which(dat$DOY > 220 & dat$Total >= 1), ]
                        # }
                        # if(species == "Leonard's Skipper"){
                        #   dat <- dat[-which(dat$DOY < 220 & dat$Total >= 1), ]
                        # }
                        
                        if(sum(dat$Total) < 20|length(unique(dat$SiteID)) < 2|length(unique(dat$Year)) < 2|
                           length(unique(dat$SiteYear))<5|length(unique(dat$RegYear))<2) {
                          mod <- NA
                        }else{
                          safe_bam <- purrr::safely(bam)
                          
                          modtime_nb <- system.time({ 
                            mod_nb <- safe_bam(Total ~
                                                 # s(zlistlength, bs = "cr") +
                                                 te(lat, lon, AccumDD, bs = c("tp", "cr"), k = c(5, 30), d = c(2, 1)) +
                                                 s(SiteYear, bs = "re") +
                                                 s(Year, bs = "re") +
                                                 s(SiteID, bs = "re") +
                                                 s(RegYear, DOY, bs = "fs", k = 5, m = 1),
                                               family = nb(theta = NULL, link = "log"),
                                               # family = poisson(link = "log"),
                                               data = dat,
                                               method = "fREML",
                                               discrete = TRUE, nthreads = 2,
                                               control = list(maxit = 500))
                            
                          })
                          
                          modtime_po <- system.time({ 
                            mod_po <- safe_bam(Total ~
                                                 # s(zlistlength, bs = "cr") +
                                                 te(lat, lon, AccumDD, bs = c("tp", "cr"), k = c(5, 30), d = c(2, 1)) +
                                                 s(SiteYear, bs = "re") +
                                                 s(Year, bs = "re") +
                                                 s(SiteID, bs = "re") +
                                                 s(RegYear, DOY, bs = "fs", k = 5, m = 1),
                                               # family = nb(theta = NULL, link = "log"),
                                               family = poisson(link = "log"),
                                               data = dat, 
                                               method = "fREML",
                                               discrete = TRUE, nthreads = 2,
                                               control = list(maxit = 500))
                            
                          })
                          
                          if(is.null(mod_po$error) == TRUE & is.null(mod_nb$error) == TRUE){
                            # compare models
                            if(AIC(mod_po) <= AIC(mod_nb)){
                              mod <- mod_po
                              modtime <- modtime_po
                            }else{
                              mod <- mod_nb
                              modtime <- modtime_nb  
                            }
                          }else{
                            mod <- ifelse(is.null(mod_po$error) == TRUE, mod_po, mod_nb)
                            modtime <- ifelse(is.null(mod_po$error) == TRUE, modtime_po, modtime_nb)
                          }
                        }
                      }
                      
                      if(is.null(mod$error)){
                        pars$modtime <- as.numeric(modtime)[1]
                        pars$AIC <- AIC(mod$result)
                        summod <- summary(mod$result)
                        pars$N <- summod$n
                        pars$dev.expl <- summod$dev.expl
                        outlist <- list()
                        outlist[["params"]] <- pars
                        outlist[["gammod"]] <- mod$result
                        outlist[["gamtime"]] <- as.numeric(modtime)[1]
                        outlist[["datGAM"]] <- dat
                        saveRDS(outlist, paste(species, "final", "rds", sep = "."))
                        
                      }else{
                        outlist <- list()
                        outlist[["params"]] <- pars
                        outlist[["gammod"]] <- mod$error
                        outlist[["datGAM"]] <- dat
                        saveRDS(outlist, paste("gamerr", species, "rds", sep = "."))
                      }
                      return(sp)
                    }

if(.Platform$OS.type == "windows"){
  stopCluster(cl)
}


mod <- mod_po
pltdat <- dat %>% filter(SiteID == "034", Year == "2000")

preddat <- gdd %>% 
  ungroup() %>% 
  mutate(Year = as.factor(as.character(Year))) %>% 
  filter(SiteID == "034", Year == "2000", DOY <= 300 & DOY >= 90) %>% 
  mutate(RegYear = paste(region, Year, sep = "_"),
         SiteYear = paste(SiteID, Year, sep = "_"),
         zlistlength = 0)

# prediction with simulated counts, stochastic but integers (if n is odd)
Xp <- predict.gam(object = mod, newdata = preddat, type="lpmatrix") ## map coefs to fitted curves
beta <- coef(mod)
Vb   <- vcov(mod) ## posterior mean and cov of coefs
n <- 100 # choose number of simulations
mrand <- MASS::mvrnorm(n, beta, Vb) ## simulate n rep coef vectors from posterior
ilink <- family(mod)$linkinv
linklppreds <- Xp %*% t(mrand)
nbpreds <- apply(X = linklppreds,
                 MARGIN = 1,
                 FUN = function(x){
                   # temp <- sort(x)
                   # bounds <- quantile(1:n, probs = c(0.025, 0.975))
                   # x <- temp[bounds[1]:bounds[2]]
                   x <- ilink(x)
                   x <- rnbinom(n = length(x),
                                mu = x,
                                size = mod$family$getTheta(TRUE))
                   # x <- quantile(x, .5)
                   return(x)
                 })
preddat$adjY <- predict.gam(object = mod, newdata = preddat, type="response") 



# PopIndex via UKBMS (missing counts imputed, trapezoid by week)



# PopIndex via GAM predictions (all weekly counts imputed, trapezoid by week, uncertainty quantified)





