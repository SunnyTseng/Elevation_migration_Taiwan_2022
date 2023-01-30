############################################################################
### Author: Sunny                                                        ###
### Project: RADI                                                        ###
### Input: file name of the downloaded eBird Basic Dataset (EBD)         ###
### Output: none, but save the claned file in the folder                 ###
############################################################################

data_preparation_ebird <- function(EBD = here("data", 
                                              "raw", 
                                              "ebd_TW_relFeb-2022", 
                                              "ebd_TW_relFeb-2022.txt"), 
                                   path = here("data", 
                                               "analyzed", 
                                               "data_eBird_qualified.csv"),
                                   start_year = 2014,
                                   end_year = 2021,
                                   effort_max_distance = 10,
                                   effort_max_duration = 300){
  
  data_filter <- EBD %>%
    auk_ebd() %>%
    auk_complete() %>%
    auk_distance(range(0, effort_max_distance)) %>%
    auk_duration((range(0, effort_max_duration))) %>%
    auk_protocol(c("Traveling", "Stationary")) %>%
    auk_year(seq(start_year, end_year)) %>%
    auk_filter(file = tempfile()) 
  
  data_select <- data_filter$output %>%
    auk_ebd() %>%
    auk_select(file = tempfile(),
               select = c("global unique identifier", "common name", "scientific name", "observation count",
                          "latitude", "longitude", "observation date", "time observations started", "observer id", 
                          "sampling event identifier", "protocol type", "protocol code", "duration minutes",
                          "effort distance km", "number observers", "all species reported",
                          "group identifier"))
  
  data <- fread(data_select)
  
  #data <- read_ebd(data_select) 
  # automatically runs auk_unique() to remove shard checklists
  # automatically runs auk_rollup() to remove subspecies
  
  data_cleaned <- data %>%
    as_tibble() %>%
    clean_names() %>%
    mutate(filt = paste0(common_name, group_identifier)) %>%
    mutate(year = observation_date %>% year(),
           month = observation_date %>% month(),
           week = observation_date %>% week(),
           day = observation_date %>% yday(),
           hour = time_observations_started %>% hms() %>% hour(),
           observation_count = observation_count %>% as.numeric()) %>%
    filter(group_identifier == "" | !duplicated(filt)) %>%
    drop_na(observation_count)

write_csv(data_cleaned, path)
return(NULL)
}
