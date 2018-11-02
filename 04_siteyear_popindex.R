source('01_data_prep.R')



tmp <- readRDS("Appalachian Eyed Brown.final.rds")
mod <- tmp$gammod
dat <- tmp$datGAM
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



# PopIndex via UKBMS (missing counts imputed, trapezoid by week)




# PopIndex via GAM predictions (all weekly counts imputed, trapezoid by week, uncertainty quantified)





