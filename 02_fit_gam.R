# ON hold until I get new data through 2016
# then rerun GAMs and do population index all at once

#script combining:
#counts filtering
#gam predictions for a collated index
#differs from poptrends.R by accounting for zero counts at sites
# use sites where sum(total) >= 5, YearSeen at site > 1
# otherwise Site's population indices very low and bring down collated index

source('01_data_prep.R')


# Fit a generalized additive model (GAM) for each species
# This model will be used to impute missing counts a la UKBMS
# OR be used to compute site x year population indices


# run a subset for testing
# allspecies <- allspecies[c(81:84), ]

ncores <- 4
if(ncores > (parallel::detectCores() / 2)){
  ncores <- parallel::detectCores() / 2
}
cl <- makePSOCKcluster(ncores)
registerDoParallel(cl)

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
                               YearTotal = sum(Total))
                      
                      # what if cutoff for inclusion is really open?
                      # could filter out sites with lower effort later
                      # also this better accounts for SiteYear zero counts, 
                      # which are important for population trends
                      dat <- counts %>% filter(YearTotal >= 0, SurvPerYear >= 5)
                      
                      mod <- list()
                      
                      if(nrow(dat) == 0){
                        mod$error <- "no data"
                      }else{
                        
                        dat$Year <- as.factor(as.character(dat$Year))
                        dat$region <- as.factor(as.character(dat$region))
                        dat$SiteID <- as.factor(as.character(dat$SiteID))
                        dat$SiteYear <- as.factor(paste(dat$SiteID, dat$Year, sep = "_"))
                        dat$zlistlength[which(is.na(dat$zlistlength))] <- 0
                        dat$RegYear <- as.factor(paste(dat$region, dat$Year, sep = "_"))
                        dat <- as.data.frame(dat)
                        
                        dat <- dat[which(!is.na(dat$AccumDD)), ]
                        
                        # #silly filters for univoltine species with outliers
                        # if(species == "Baltimore"){
                        #   dat <- dat[-which(dat$DOY > 220 & dat$Total >= 1), ]
                        # }
                        # if(species == "Leonard's Skipper"){
                        #   dat <- dat[-which(dat$DOY < 220 & dat$Total >= 1), ]
                        # }
                        
                        if(sum(dat$Total) < 20|length(unique(dat$SiteID)) < 2|length(unique(dat$Year)) < 2|
                           length(unique(dat$SiteYear))<5|length(unique(dat$RegYear))<2) {
                          mod$error <- "few data"
                        }else{
                          safe_bam <- purrr::safely(bam)
                          
                          modtime_nb <- system.time({ 
                            mod_nb <- safe_bam(Total ~
                                                 # s(zlistlength, bs = "cr") +
                                                 te(lat, lon, AccumDD, bs = c("tp", "cc"), k = c(5, 30), d = c(2, 1)) +
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
                            if(AIC(mod_po$result) <= AIC(mod_nb$result)){
                              mod <- mod_po
                              modtime <- modtime_po
                            }else{
                              mod <- mod_nb
                              modtime <- modtime_nb  
                            }
                          }else{
                            if(is.null(mod_po$error) == TRUE){
                              mod <- mod_po
                              modtime <- modtime_po
                            }else{
                              mod <- mod_nb
                              modtime <- modtime_nb
                            }
                          }
                        }
                      }
                      
                      if(is.null(mod$error)){
                        pars$modtime <- as.numeric(modtime)[1]
                        pars$AIC <- AIC(mod$result)
                        summod <- summary(mod$result)
                        pars$N <- summod$n
                        pars$dev.expl <- summod$dev.expl
                        pars$family <- summod$family$family
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

