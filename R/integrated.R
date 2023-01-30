#####################
### Author: Sunny ###
### Project: RADI ###
#####################

###############
### Library ###
###############
# data management
library(tidyverse) #
library(data.table) #
library(here) #
library(auk) #
library(janitor) #
library(furrr) #
library(ranger) #
library(Boruta) #
library(mgcv) #

# GIS related 
library(raster) #
library(rgdal) #  
library(sf) #
library(twmap) #

# plot related
library(RColorBrewer) #
library(lattice) #
library(ggcorrplot) #
library(plotly) #
library(cowplot) #
library(PupillometryR) #

#################
### Functions ###
#################
# resolve namespace conflicts
select <- dplyr::select
map <- purrr::map
projection <- raster::projection

# functions for data preparation
source(here("R", "data_preparation_ebird.R"))
source(here("R", "data_preparation_predictors.R"))
source(here("R", "data_preparation_prediction_surface.R"))
source(here("R", "data_preparation_target_species.R"))

# functions for modelling
source(here("R", "modelling_random_forest.R"))
source(here("R", "modelling_stixel_grouping.R"))
source(here("R", "modelling_GAM.R"))
source(here("R", "modelling_evaluation.R"))
source(here("R", "prediction_maps.R"))

########################
### Data Preparation ###
########################
### cleaned data for eBird and predictors, umcommend when needed
data_preparation_ebird(EBD = here("data", 
                                  "raw", 
                                  "ebd_TW_relFeb-2022", 
                                  "ebd_TW_relFeb-2022.txt"), 
                       path = here("data", 
                                   "analyzed", 
                                   "data_eBird_qualified.csv"),
                       start_year = 2014,
                       end_year = 2021,
                       effort_max_distance = 10,
                       effort_max_duration = 300)

### cleaned predictors data according to the eBird checklists, uncommend when needed
data_preparation_predictors(dir_tiff = here("data",
                                            "raw" , 
                                            "Taiwan_environmental_dataset-master", 
                                            "GeoTIFF_unzip"),
                            dir_eBird = here("data", 
                                             "analyzed", 
                                             "data_eBird_qualified.csv"),
                            path = here("data", 
                                        "analyzed", 
                                        "data_eBird_qualified_predictors.csv"))

### prediction surface for making predictions, uncommend when needed
data_preparation_prediction_surface(dir_tiff = here("data", 
                                                    "raw", 
                                                    "Taiwan_environmental_dataset-master", 
                                                    "GeoTIFF_unzip"),
                                    path_data_frame = here("data", 
                                                           "analyzed", 
                                                           "prediciton_surface.csv"),
                                    path_tif = here("data", 
                                                    "analyzed", 
                                                    "prediction_surface.tif"))

### define target species and save the data for target species, uncommend when needed
data_preparation_target_species(dir_eBird = here("data", 
                                                 "analyzed", 
                                                 "data_eBird_qualified.csv"),
                                dir_predictors = here("data", 
                                                      "analyzed", 
                                                      "data_eBird_qualified_predictors.csv"),
                                target_species = "Yuhina brunneiceps",
                                path = here("data", 
                                            "analyzed", 
                                            paste0("data_eBird_qualified_combined_", target_species, ".csv")))


##########################
### Variable selection ###
##########################

target_species <- "Yuhina brunneiceps"

# import data

data <- read_csv(here("data", 
                      "analyzed", 
                      paste0("data_eBird_qualified_combined_", target_species, ".csv")))

data <- data %>%
  mutate(detection = detection %>% as.factor(),
         protocol_type = protocol_type %>% as.factor(),
         other_proad = other_proad %>% as.factor(),
         year = year %>% as.factor()) %>%
  mutate(effort_distance_km = if_else(is.na(effort_distance_km), 0, effort_distance_km)) %>%
  drop_na()

# sub sampling

data_sub <- data %>%
  # read in data and save the lat lon info
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  mutate(longitude = unlist(map(.$geometry,1)),
         latitude = unlist(map(.$geometry,2))) %>%
  # filter the data points collected outside of main island
  st_transform(3824) %>%
  st_join(tw_county, join = st_within) %>%
  filter(!COUNTYID %in% c("Z", "W", "X", NA)) %>%
  # create grid and perform sub-sampling
  st_transform(3826) %>%
  mutate(y = unlist(map(.$geometry,1)),
         x = unlist(map(.$geometry,2))) %>% 
  mutate(y_cut = cut(y, breaks = seq(from = (range(.$y)[1] - 2000),
                                     to = (range(.$y)[2] + 2000),
                                     by = 2000)),
         x_cut = cut(x, breaks = seq(from = (range(.$x)[1] - 2000),
                                     to = (range(.$x)[2] + 2000),
                                     by = 2000))) %>%
  unite(cell, y_cut:x_cut, remove = FALSE) %>%
  group_by(detection, week, cell) %>%
  sample_n(size = 1) %>%
  ungroup() %>%
  as_tibble() %>%
  select(-y_cut, -x_cut, -geometry)

count(data_sub, detection) %>% 
  mutate(percent = n / sum(n))

distribution <- tw_county %>%
  st_crop(xmin = 119, xmax = 123, ymin = 20, ymax = 26) %>%
  ggplot() +
  geom_sf(size = 0.2) +
  geom_point(data = filter(data_sub, detection == 0), aes(x = longitude, y = latitude), size = 0.05, colour = "grey39", alpha = 0.4) +
  geom_point(data = filter(data_sub, detection == 1), aes(x = longitude, y = latitude), size = 0.05, colour = "forestgreen", alpha = 0.4) +
  theme_bw() 

elevation_raw <- data_sub %>%
  filter(detection == 1) %>%
  group_by(week) %>%
  summarise(elemean = dtm_ele %>% mean(),
            elesd = dtm_ele %>% sd(),
            checklists = n()) %>%
  mutate(elese = elesd/sqrt(checklists)) %>%
  ggplot(aes(x=week, y=elemean)) + 
  geom_line(col = "#339966", size = 0.8) +
  geom_point(size = 2) +
  geom_area(fill = "#339966", alpha = 0.5) +
  geom_errorbar(aes(ymin = elemean - elese, ymax = elemean + elese), width=.2,
                position=position_dodge(0.05)) +
  theme_bw()



# predictors selection, would take a while :)
predictors <- modelling_random_forest(data = data_sub, 
                                      method = "ranger",
                                      cor_threshold = 0.8,
                                      max_vars = 20)

#####################################
### Data split and create stixels ###
#####################################
stixels <- modelling_stixel_grouping(data = data_sub,
                                     predictors = predictors,
                                     split = 0.8,
                                     temporal_resolution = 7,
                                     stixel_height = 40)

##################################
### Model training and testing ###
##################################
set.seed(100)
nb <- modelling_GAM(stixels = stixels, family = "nb", predictors = predictors, workers = 16)
ziplss <- modelling_GAM(stixels = stixels, family = "ziplss", predictors = predictors, workers = 16)
models <- inner_join(nb, ziplss %>% select(day_of_year, m_ziplss), by = "day_of_year")
rm(nb)
rm(ziplss) # Time difference of 29.26519 mins

models_test <- modelling_evaluation(models = nb, family = "nb", workers = 16) 


#######################################
### Prediction using selected model ###
### Visualization of maps           ###
#######################################

models_map_nb <- prediction_maps(models = nb,
                                 family = "nb",
                                 workers = 16,
                                 duration_minutes = 60,
                                 effort_distance_km = 1, 
                                 number_observers = 1, 
                                 hour = 7,
                                 dir_pred_surf = here("data", "analyzed", "prediciton_surface.csv"),
                                 dir_pred_tif = here("data", "analyzed", "prediction_surface.tif"),
                                 quantile = 0.2)
                                     
### Abundance maps
for(i in 1:53){
  png(here("docs", target_species,
           paste0(target_species, "_week_", i, ".png")), 
      res = 300, width = 3, height = 4, units = 'in')
  print(models_map_nb$map_nb[[i]])
  dev.off()
}

