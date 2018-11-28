source('01_data_prep.R')

modfiles <- list.files(path = "gams", full.names = TRUE)
modfiles <- modfiles[grep(pattern = "final", x = modfiles, fixed = TRUE)]
sp_params <- list()

for (mf in seq_along(modfiles)){
  
  tmp <- readRDS(modfiles[mf])
  mod <- tmp$gammod
  counts <- tmp$datGAM
  sp_params[[mf]] <- tmp$params

    
  latin <- allspecies %>% 
    ungroup() %>% 
    filter(CommonName == tmp$params$CommonName) %>% 
    select(CombinedLatin) %>% 
    as.character()
  
  # plot GAM predictions for species/model
  preds <- gdd %>% 
    mutate(SiteYear = paste(SiteID, Year, sep = "_")) %>% 
    filter(SiteYear %in% unique(counts$SiteYear)) %>%
    group_by(SiteID, Year) %>% 
    mutate(SiteYearGDD = max(AccumDD)) %>% 
    filter(DOY %in% seq(90, 305, 1)) %>% 
    ungroup() %>% 
    mutate(RegYear = as.factor(paste(region, Year, sep = "_" )),
           Year = as.factor(as.character(Year)),
           SiteID = as.factor(SiteID),
           SiteYear = as.factor(SiteYear))
  
  preds$adjY <- predict.gam(object = mod, newdata = preds, type="response")
  
  preds <- preds %>% 
    group_by(SiteYear) %>%
    mutate(Gamma = adjY / as.vector(sum(adjY)),
           SiteYearTotal = sum(adjY)) %>%
    ungroup() %>% 
    filter(adjY > 0) %>% 
    # mutate(ztotal = ((.7 - .1) * (SiteYearTotal - min(SiteYearTotal))
    #                  /(max(SiteYearTotal) - min(SiteYearTotal))) + .1) %>% 
    group_by(region) %>% 
    filter(SiteYear %in% sample(unique(SiteYear), 10, replace = TRUE))
  
  preds$region <- factor(preds$region, levels = c("NW", "NE", "SW", "CN"))
  
  # outliers
  outs <- counts %>% 
    filter(Total > 0) %>% 
    dplyr::select(AccumDD, DOY, region) %>% 
    filter(complete.cases(.))
  
  outs$region <- factor(outs$region, levels = c("NW", "NE", "SW", "CN"))
  
  # # if(tmp$params$model == "doy"){
  # gamplt <- ggplot(preds, aes(x = DOY, y = Gamma, group = SiteYear, color = SiteYearGDD)) +
  #   geom_path(alpha = .5) + 
  #   scale_color_viridis() + 
  #   facet_wrap(~region, ncol = 2) +
  #   geom_rug(data = outs, aes(x = DOY), sides="b", inherit.aes = FALSE, alpha = 0.3) +
  #   ggtitle(paste0(tmp$params$CommonName, " (", latin, ")"),
  #           subtitle = "seasonal phenology modeled on calendar scale") +
  #   labs(color = "Total degree-days\n for site and year") +
  #   labs(x = "Day of year") +
  #   labs(y = "Scaled phenology (model predictions)")
  # ggsave(filename = paste(tmp$params$CommonName, "GAM", "DOY", "png", sep = "."), 
  #        plot = gamplt, device = "png", path = "plots", width = 8, height = 6, units = "in")
  # # }else{
  gamplt <- ggplot(preds, aes(x = AccumDD, y = Gamma, group = SiteYear, color = SiteYearGDD)) +
    geom_path(aes(alpha = log(SiteYearTotal))) + 
    scale_color_viridis() + 
    facet_wrap(~region, ncol = 2) +
    geom_rug(data = outs, aes(x = AccumDD), sides="b", inherit.aes = FALSE,  alpha = 0.3) +
    ggtitle(paste0(tmp$params$CommonName, " (", latin, ")"),
            subtitle = "seasonal phenology modeled on degree-day scale") +
    labs(color = "Total degree-days\n for site and year") +
    labs(x = "Degree-days accumulated (10/30C thresholds)") +
    labs(y = "Scaled phenology (model predictions)")
  # }
  ggsave(filename = paste(tmp$params$CommonName, "GAM", "GDD", "png", sep = "."), 
         plot = gamplt, device = "png", path = "plots", width = 8, height = 6, units = "in")
  
  
}

outparams <- bind_rows(sp_params)
saveRDS(outparams, "modparams.5.rds")




modfiles <- list.files()
modfiles <- modfiles[grep(pattern = "final", x = modfiles, fixed = TRUE)]
sp_params <- list()

for (mf in seq_along(modfiles)){
  
  tmp <- readRDS(modfiles[mf])
  sp_params[[mf]] <- tmp$params
}
outparams <- bind_rows(sp_params)

outparam10 <- readRDS("modparams.rds")

