source('01_data_prep.R')

modfiles <- list.files("gams", full.names = TRUE)
modfiles <- modfiles[grep(pattern = "final", x = modfiles, fixed = TRUE)]
sp_pops <- list()

for (mf in seq_along(modfiles)){
  # for (mf in 1:51){
    
  tmp <- readRDS(modfiles[mf])
  mod <- tmp$gammod
  dat <- tmp$datGAM
  
  
  
  # PopIndex via UKBMS (missing counts imputed, trapezoid by week)
  # UKBMS count index approximation
  TrapezoidIndex <- function(timescale, counts){
    dat <- data.frame(t = timescale, y = counts)
    dat <- arrange(dat, timescale)
    if(nrow(dat) == 1){
      return(dat$y)
    }else{
      temp <- rep(NA, nrow(dat)-1)
      for (i in 2:length(dat$t)){
        temp[i-1] <- (dat$y[i] + dat$y[i-1]) * (dat$t[i] - dat$t[i-1]) / 2
      }
    }
    return(sum(temp))
  }
  
  surv <- dat %>% 
    dplyr::select(SiteID, Week, Year, SiteYear) %>% 
    distinct()
  allsurv <- surv %>% 
    complete(SiteID, Week, Year) %>% 
    mutate(SiteYear = paste(SiteID, Year, sep = "_")) %>% 
    filter(Week <= 30) %>% 
    filter(SiteYear %in% unique(surv$SiteYear))
  missing <- anti_join(allsurv, surv)
  
  impute_days <- expand.grid(Year = unique(missing$Year), Week = c(1:30)) %>% 
    mutate(DOY = ifelse(as.numeric(as.character(Year)) %% 4 == 0, 92 + Week * 7 - 4, 91 + Week * 7 - 4)) # midweek to impute
  
  missing <- missing %>% 
    left_join(impute_days) %>% 
    mutate(SiteDate = as.Date(paste(DOY, Year, sep = "-"), format = "%j-%Y")) %>% 
    left_join(gdd[, c("SiteID", "SiteDate", "lat", "lon", "AccumDD", "region")], by = c("SiteID", "SiteDate")) %>% 
    mutate(RegYear = paste(region, Year, sep = "_")) %>% 
    mutate(Total = predict(mod, newdata = ., type = "response"),
           datatype = "imputed")
  
  impute_all <-  expand.grid(Year = unique(allsurv$Year), Week = c(1:30)) %>% 
    mutate(DOY = ifelse(as.numeric(as.character(Year)) %% 4 == 0, 92 + Week * 7 - 4, 91 + Week * 7 - 4)) # midweek to impute
  
  allpred <- allsurv %>% 
    left_join(impute_all) %>% 
    mutate(SiteDate = as.Date(paste(DOY, Year, sep = "-"), format = "%j-%Y")) %>% 
    left_join(gdd[, c("SiteID", "SiteDate", "lat", "lon", "AccumDD", "region")], by = c("SiteID", "SiteDate")) %>% 
    mutate(RegYear = paste(region, Year, sep = "_")) %>% 
    mutate(Total = predict(mod, newdata = ., type = "response"),
           datatype = "predicted")
  
  # check for boundary problems
  species <- tmp$params$CommonName
  jpeg(paste0(species, ".jpg"), width = 1000, height = 1000, quality = 75)
  par(mfcol=c(3, 2), cex = 1)
  plot(missing$AccumDD, missing$Total, main = paste("imputed", species, sep = "_"))
  plot(dat$AccumDD, dat$Total, main = paste("observed", species, sep = "_"))
  plot(allpred$AccumDD, allpred$Total, main = paste("predicted", species, sep = "_"))
  plot(missing$DOY, missing$Total, main = paste("imputed", species, sep = "_"))
  plot(dat$DOY, dat$Total, main = paste("observed", species, sep = "_"))
  plot(allpred$DOY, allpred$Total, main = paste("predicted", species, sep = "_"))
  dev.off()

  alldat <- dat %>%
    mutate(datatype = "observed") %>%
    bind_rows(missing) %>%
    group_by(SiteYear, SiteID, Year, lat, lon, RegYear) %>%
    summarise(Index = TrapezoidIndex(DOY, Total)) %>%
    mutate(CommonName = species,
           method = "ukbms")

  covs <- dat %>% 
    dplyr::select(SiteYear, Total, listlength, duration, SurvPerYear, YearTotal) %>% 
    group_by(SiteYear) %>% 
    summarise(meanLL = mean(listlength, na.rm = TRUE),
              meanTime = mean(duration, na.rm = TRUE),
              SurvPerYear = SurvPerYear[1],
              YearTotal = YearTotal[1],
              posWeeks = length(which(Total > 0)))
  
  alldat <- alldat %>% left_join(covs)
  
  allpreds <- allpred %>%
    group_by(SiteYear, SiteID, Year, lat, lon, RegYear) %>%
    summarise(Index = TrapezoidIndex(DOY, Total)) %>%
    mutate(CommonName = species,
           method = "gampred")
  
  allpreds <- allpreds %>% left_join(covs)
  
  out <- bind_rows(alldat, allpreds)
  
  sp_pops[[mf]] <- out
}

allpops <- bind_rows(sp_pops)
saveRDS(allpops, "allpops.5.rds")

# NOTES ON PLOTS CHECKING PREDICTIONS/IMPUTED MISSING COUNTS

# Checkered White has way too high values from poor model imputation
# Others notes from plots
# Red-banded Hairstreak imputed values in 2nd generation higher than observed
# Northern Metalmark has imputed 10x maximum observed
# Long Dash, Little Wood Satyr, Hobomok Skipper, Falcate Orangetip, Clouded Sulphur, 
# have imputed 1.5x maximum observed


# TODO: Uncertainty about GAM predictions
# PopIndex via GAM predictions (all weekly counts imputed, trapezoid by week, uncertainty quantified)


# example with one SiteYear
pltdat <- dat %>% filter(SiteID == "097", Year == "2007")

preddat <- gdd %>% 
  ungroup() %>% 
  mutate(Year = as.factor(as.character(Year))) %>% 
  filter(SiteID %in% unique(pltdat$SiteID), Year %in% unique(pltdat$Year), DOY <= 300 & DOY >= 90) %>% 
  mutate(RegYear = paste(region, Year, sep = "_"),
         SiteYear = paste(SiteID, Year, sep = "_"),
         zlistlength = 0)

# prediction with simulated counts, stochastic but integers (if n is odd)
Xp <- predict.gam(object = mod, newdata = preddat, type="lpmatrix") ## map coefs to fitted curves
beta <- coef(mod)
Vb   <- vcov(mod) ## posterior mean and cov of coefs
n <- 500 # choose number of simulations
mrand <- MASS::mvrnorm(n, beta, Vb) ## simulate n rep coef vectors from posterior
ilink <- family(mod)$linkinv
linklppreds <- Xp %*% t(mrand)
# nbpreds <- apply(X = linklppreds,
#                  MARGIN = 1,
#                  FUN = function(x){
#                    # temp <- sort(x)
#                    # bounds <- quantile(1:n, probs = c(0.025, 0.975))
#                    # x <- temp[bounds[1]:bounds[2]]
#                    x <- ilink(x)
#                    x <- rnbinom(n = length(x),
#                                 mu = x,
#                                 size = mod$family$getTheta(TRUE))
#                    # x <- quantile(x, .5)
#                    return(x)
#                  })
nbpreds <- ilink(linklppreds)
matplot(nbpreds, type = "l")


preddat$adjY <- predict.gam(object = mod, newdata = preddat, type="response") 





