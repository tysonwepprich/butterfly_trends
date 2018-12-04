# take abundance trends and see what site factors affect mean species trajectory
source("01_data_prep.R")
# Statewide pesticide use by crop category
pest <- read.delim("data/HighEstimate_AgPestUsebyCropGroup92to16.txt", header = TRUE, sep = "\t") %>% 
  filter(State == "Ohio") %>% 
  tidyr::gather(Crop, KG, Corn:Other_crops) %>% 
  # filter(Crop == "Soybeans") %>% 
  group_by(Compound) %>% 
  mutate(Total = sum(KG, na.rm = TRUE)) %>% 
  # filter(Total > 10000) %>% 
  droplevels()

plt <- ggplot(pest, aes(x = Year, y = KG, group = Crop)) +
  geom_point() + 
  geom_smooth(method = "glm", se = TRUE, method.args = list(family = "poisson")) +
  facet_wrap(~Compound, scales = "free_y") +
  theme_bw()
plt


# this is old MOA/pesticide traits file, might need updating
download.file("https://raw.githubusercontent.com/slomascolo/NCEAS-RENCI_2014/master/Pesticides/PestGrouping/OHpest_traits_allRACs.csv", 
              destfile = "data/moa.csv", method = "curl")
moa <- read.csv("data/moa.csv", header = TRUE)

pest_moa <- pest %>% 
  dplyr::select(Compound, Total) %>% 
  distinct() %>% 
  left_join(moa, by = c("Compound" = "COMPOUND"))
write.csv(pest_moa, "data/moa_update.csv", row.names = FALSE)



library(raster)
r <- raster("data/CDL/CDL_2016_clip_20181130155428_520910358.tif")


