#################################
### Library and .RData import ###
#################################
library(tidyverse)
library(purrr)
library(here)

load(here("WESI.RData"))

###
### objective 0: model evaluation - performance seasonality
###

# evaluation <- models_test %>%
#   transmute(cor_spearman = purrr::map_dbl(.x = test_pred_nb,
#                                           .f =~ cor(.x$abd, .x$observation_count, method = "spearman")),
#             cor_pearson = purrr::map_dbl(.x = test_pred_nb,
#                                           .f =~ cor(.x$abd, .x$observation_count, method = "pearson")),
#             cor_kendall = purrr::map_dbl(.x = test_pred_nb,
#                                          .f =~ cor(.x$abd, .x$observation_count, method = "kendall")),
#             test_size = purrr::map_dbl(.x = test_pred_nb, .f =~ nrow(.x))) %>%
#   mutate(day_of_year = seq(from = 1, by = 7, length.out = nrow(models_test)),
#          species = "WTRO") 
#write_csv(evaluation, here("data", "analyzed", "model_evaluation_Myiomela leucura.csv"))

eval_WESI <- read_csv(here("data", "analyzed", "model_evaluation_Heterophasia auricularis.csv"))
eval_WTRO <- read_csv(here("data", "analyzed", "model_evaluation_Myiomela leucura.csv"))
eval_TAYU <- read_csv(here("data", "analyzed", "model_evaluation_Yuhina brunneiceps.csv"))
eval_all <- eval_WESI %>%
  rbind(eval_WTRO) %>%
  rbind(eval_TAYU)

evaluation_fig <- eval_all %>%
  ggplot(aes(colour = species, group = species)) +
  geom_line(aes(x = day_of_year, y = cor_spearman), size = 1.1) +
  geom_point(aes(x = day_of_year, y = cor_spearman), size = 2, shape = 16) +
  geom_line(aes(x = day_of_year, y = (test_size - min(test_size))/(max(test_size) - min(test_size))),
            size = 1.1, linetype = 2, alpha = 0.6) +
  scale_color_brewer(palette = "Dark2") +
  scale_y_continuous(
    # Features of the first axis
    name = "Spearman's correlation",
    # Add a second axis and specify its features
    sec.axis = sec_axis(trans =~.*784 + 1173 ,
                        name = "# of test observations")) +
  labs(x = "Day of year") +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 12),
        axis.title.y.left = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.y.right = element_text(margin = margin(t = 0, r = 0, b = 0, l = 10)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) 

evaluation_fig


### Radar chart
test <- evaluation_WESI %>%
  mutate("1" = 0.6, "2" = 0) %>%
  select("1", "2", cor_pearson, cor_kendall, cor_spearman) %>%
  t() %>%
  as.data.frame()
colnames(test) <- evaluation_data$day_of_year
rownames(test) <- c("1", "2", "cor_pearson", "cor_kendall", "cor_spearman")

radarchart(test)


###
### objective 0: model evaluation - coefficient table
###

# ggplot function
plot_gam <- function(m, title = NULL, ziplss = c("presence", "abundance")) {
  # capture plot
  tmp <- tempfile()
  png(tmp)
  p <- plot(m, pages = 1)
  dev.off()
  unlink(tmp)
  
  # drop addition models in ziplss
  if (m$family$family == "ziplss") {
    is_presence <- map_lgl(p, ~ str_detect(.$ylab, "^s\\.1"))
    if (ziplss == "presence") {
      p <- p[is_presence]  
    } else {
      p <- p[!is_presence]
    }
  }
  
  # extract data
  p_df <- map_df(p, ~ tibble(cov = rep(.$xlab, length(.$x)),
                             x = .$x, fit = .$fit, se = .$se))
  
  # plot
  g <- ggplot(p_df) +
    aes(x = x, y = fit,
        ymin = fit - se, ymax = fit + se) +
    geom_ribbon(fill = "grey80") +
    geom_line(col = "blue") +
    facet_wrap(~ cov, scales = "free_x") +
    labs(x = NULL,
         y = "Smooth function",
         title = title)
  print(g)
  invisible(p_df)
}

plot_gam(nb$m_nb[[5]], title = "Negative Binomial GAM")


###
### Objective 1: fine scale information - weekly change of population mean elevation, area coverage
###
dtm <- here("data", 
            "raw", 
            "Taiwan_environmental_dataset-master", 
            "GeoTIFF_unzip", 
            "GeoTIFF_epsg3826_DTM", 
            "G1km_TWD97-121_DTM_ELE.tif") %>%
  raster()

area_change <- models_map_nb %>%
  mutate(area = map2_dbl(.x = threshold, .y = map_pred_nb, .f = ~ .y %>% 
                           filter(abd > .x) %>%
                           dim() %>%
                           .[1]),
         area_diff = (area- lag(area))/area) %>%
  mutate(season = cut(day_of_year, 
                      breaks = c(1 , 64 , 106, 253, 344, Inf),
                      right = FALSE,
                      labels = c("NB", "BB", "BR", "AB", "NB"))) 

annual_ele <- area_change %>%
  mutate(ele_data = map2(.x = map_pred_nb, .y = threshold, .f = ~ .x %>%
                           mutate(abd = if_else(abd > .y, abd, 0),
                                  ele = extract(dtm, .x %>% select(x, y) %>% SpatialPoints(proj4string = CRS("+init=epsg:3826"))),
                                  ele_category = cut(ele, breaks = c(seq(0, 2500, by = 250), Inf), include.lowest = TRUE),
                                  abd = ceiling(abd))),
         ele_hist = map(.x = ele_data, .f = ~ .x %>%
                          map2(.x = .$abd, .y = .$ele, .f = ~ rep(.y, times = .x)) %>% 
                          flatten_dbl() %>%
                          as_tibble(.)),
         ele_average = map_dbl(.x = ele_hist, .f = ~ .x %>% pull(value) %>% mean()),
         ele_Q1 = map_dbl(.x = ele_hist, .f = ~ .x %>% pull(value) %>% quantile(0.25)),
         ele_Q3 = map_dbl(.x = ele_hist, .f = ~ .x %>% pull(value) %>% quantile(0.75))) %>%
  select(day_of_year, season, ele_average, ele_Q1, ele_Q3) 

# area change during the year
ggplot(area_change) +
  #geom_smooth(aes(x = day_of_year, y = area), se = FALSE, span = 0.3, color = "black", linetype = "dashed") +
  geom_point(aes(x = day_of_year, y = area, fill = season), size = 5, shape = 21) +
  scale_fill_manual(values = c(NB = "#154e70", BB = "#fff2cc", BR = "#9a3c2e", AB = "#fff2cc")) +
  theme_bw()

# elevation change during the year
ggplot(annual_ele, aes(x = day_of_year, y = ele_average)) +
  #geom_smooth(aes(x = day_of_year, y = area), se = FALSE, span = 0.3, color = "black", linetype = "dashed") +
  #geom_errorbar(aes(ymin = ele_Q1, ymax = ele_Q3), width=.2) +
  geom_point(aes(fill = season), size = 5, shape = 21) +
  scale_fill_manual(values = c(NB = "#154e70", BB = "#fff2cc", BR = "#9a3c2e", AB = "#fff2cc")) +
  theme_bw()



###
### Objective 2: broad scale information - elevation distribution and coverage in breeding and non-breeding season
###

# Taiwan map
abd_season <- area_change %>%
  group_nest(season) %>%
  mutate(abd_season = map(.x = data, .f = ~ .x$map_pred_nb %>% 
                            reduce(., inner_join, by = c("x", "y")) %>%
                            mutate(average = rowMeans(select(., starts_with("abd")))) %>%
                            select(x, y, average)),
         threshold = map_dbl(.x = data, .f = ~ .x$threshold %>% mean())) 

abd_season_map <- abd_season %>%
  filter(season == "BB" | season == "NB") %>%
  select(abd_season) %>%
  map_df(.x = ., .f = ~ .x %>% reduce(., inner_join, by = c("x", "y"))) %>%
  rename(BB_red = average.y, NB_green = average.x) %>%
  mutate(BB_red = if_else(BB_red < as.numeric(abd_season[2, 4]), 0, BB_red),
         NB_green = if_else(NB_green < as.numeric(abd_season[1, 4]), 0, NB_green)) %>%
  pivot_longer(cols = c("NB_green", "BB_red"), names_to = "season") 

# Taiwan map for non-breeding season
ggplot() +
  geom_tile(data = abd_season_map, aes(x = x, y = y), fill = "#e6e6e6") +
  geom_tile(data = abd_season_map %>% filter(value != 0 & season == "NB_green"), 
            aes(x = x, y = y, fill = value)) +
  scale_fill_gradient2(trans = "log",
                       low = "#abbeca",
                       mid = "#3e7da2",
                       high = "#154e70",
                       breaks = c(0, 2, 50),
                       limits = c(from = 0.01, to = 50)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        axis.title = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank())

# Taiwan map for breeding season
ggplot() +
  geom_tile(data = abd_season_map, aes(x = x, y = y), fill = "#e6e6e6") +
  geom_tile(data = abd_season_map %>% filter(value != 0 & season == "BB_red"), 
            aes(x = x, y = y, fill = value)) +
  scale_fill_gradient2(trans = "log",
                       low = "#d8b8b4",
                       mid = "#cf695a",
                       high = "#9a3c2e",
                       breaks = c(0, 2, 50),
                       limits = c(from = 0.01, to = 50)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        axis.title = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank())

# elevation histogram
abd_season_distribution <- abd_season %>%
  mutate(ele = map(.x = abd_season, .f = ~ .x %>%
                     mutate(ele = extract(dtm, .x %>% select(x, y) %>% SpatialPoints(proj4string = CRS("+init=epsg:3826")))) %>% 
                     mutate(ele_category = cut(ele, breaks = c(seq(0, 2500, by = 250), Inf), include.lowest = TRUE)) %>%
                     group_by(ele_category) %>%
                     summarise(abd_sum = sum(average)) %>%
                     mutate(abd_sum_relative = abd_sum/sum(abd_sum))))

# elevation histogram for non-breeding season
ggplot(data = abd_season_distribution$ele[[1]], aes(x = ele_category, y = abd_sum_relative)) +
  geom_bar(stat = "identity", colour = "black", fill = "#154e70") +
  ylim(0, 0.3) +
  theme_bw()

# elevation histogram for breeding season
ggplot(data = abd_season_distribution$ele[[3]], aes(x = ele_category, y = abd_sum_relative)) +
  geom_bar(stat = "identity", colour = "black", fill = "#9a3c2e") +
  ylim(0, 0.3) +
  theme_bw()

# elevation distribution 
season_ele <- abd_season %>%
  mutate(ele_data = map2(.x = abd_season, .y = threshold, .f = ~ .x %>%
                           mutate(abd = if_else(average > .y, average, 0),
                                  ele = extract(dtm, .x %>% select(x, y) %>% SpatialPoints(proj4string = CRS("+init=epsg:3826"))),
                                  ele_category = cut(ele, breaks = c(seq(0, 2500, by = 250), Inf), include.lowest = TRUE),
                                  abd = ceiling(abd)))) %>%
  mutate(ele_hist = map(.x = ele_data, .f = ~ .x %>%
                          map2(.x = .$abd, .y = .$ele, .f = ~ rep(.y, times = .x)) %>% 
                          flatten_dbl() %>%
                          as_tibble(.))) %>%
  select(season, ele_hist) %>%
  filter(season == "BR" | season == "NB") %>%
  unnest(ele_hist) %>%
  mutate(season_dbl = if_else(season == "NB", 0.7, 0.85),
         category = "category")


# elevation distribution of breeding and non-breeding season
ggplot(data = season_ele, aes(x = season, y = value, fill = season, colour = season)) +
  #geom_flat_violin(position = position_nudge(x = .25, y = 0), adjust = 2, trim = FALSE) +
  geom_point(data = season_ele %>% sample_frac(size = 0.3),
             aes(x = category, y = value, fill = season, colour = season),
             position = position_jitter(width = .1), size = .05, shape = 16, alpha = 0.2) +
  geom_boxplot(aes(x = season_dbl + 0.25, y = value), 
               outlier.shape = NA, width = .1, colour = "BLACK", size = 0.8,
               position = position_nudge(x = .25, y = 0)) +
  ylab("Elevation") +
  xlab("Season") +
  theme_cowplot() +
  guides(fill = FALSE, colour = FALSE) +
  scale_colour_manual(values = c("#154e70", "#9a3c2e")) +
  scale_fill_manual(values = c("#154e70", "#9a3c2e")) 

# elevation distribution of breeding and non-breeding season
ggplot(data = season_ele, aes(x = category, y = value, fill = season, group = season)) +
  geom_flat_violin(aes(fill = season), position = position_nudge(x = .1, y = 0), adjust = 1.5, trim = FALSE, alpha = .5, colour = NA)+
  scale_colour_manual(values = c("#154e70", "#9a3c2e")) +
  scale_fill_manual(values = c("#154e70", "#9a3c2e")) +
  ggtitle("Figure 10: Repeated Measures Factorial Rainclouds") +
  theme_cowplot()


#####################
### Visualization ###
#####################

### distribution of checklists before and after subsampling
sub_before <- tw_county %>%
  st_crop(xmin = 119, xmax = 123, ymin = 20, ymax = 26) %>%
  ggplot() +
  geom_sf(size = 0.1) +
  geom_point(data = filter(data, detection == 0), aes(x = longitude, y = latitude), size = 0.05, colour = "grey39", alpha = 0.4) +
  geom_point(data = filter(data, detection == 1), aes(x = longitude, y = latitude), size = 0.05, colour = "forestgreen", alpha = 0.4) +
  theme_bw() 

### stixels visualization

structure <- stixels %>%
  mutate(n_training = map_dbl(.x = train_data, .f = ~ dim(.x)[1]),
         n_training_detection = map_dbl(.x = train_data, .f = ~ .x %>% filter(observation_count != 0) %>% dim() %>% .[1]),
         n_training_non_detection = map_dbl(.x = train_data, .f = ~ .x %>% filter(observation_count == 0) %>% dim() %>% .[1]),
         rate = n_training_detection/n_training)


### model evaluation
model_eval <- models_test %>%
  mutate(RMSE_train = map_dbl(.x = models_test$m_nb, .f = ~ .x$residuals %>% .**2 %>% mean() %>% sqrt()),
         RMSE_test = map_dbl(.x = models_test$test_pred_nb, .f = ~ (.x$abd - .x$observation_count) ** 2 %>% mean() %>% sqrt()),
         deviance_explained = map_dbl(.x = m_nb, .f = ~ (.x$deviance/.x$null.deviance)),
         n_training = map_dbl(.x = train_data, .f = ~ dim(.x)[1]),
         n_training_detection = map_dbl(.x = train_data, .f = ~ .x %>% filter(observation_count != 0) %>% dim() %>% .[1]),
         n_training_non_detection = map_dbl(.x = train_data, .f = ~ .x %>% filter(observation_count == 0) %>% dim() %>% .[1]))

ggplot(data = model_eval, aes(x = day_of_year)) +
  geom_line(aes(y = n_training_non_detection), size = 1.5, colour = "steelblue4") +
  geom_line(aes(y = n_training_detection*20), size = 1.5, colour = "steelblue1") +
  #geom_line(aes(y = RMSE_test*10000, colour = "Testing"), size = 1, colour = "red", linetype = "dashed") +
  scale_y_continuous(sec.axis = sec_axis(~./20, name = "Detection")) +
  theme_bw()

ggplot(data = model_eval, aes(x = day_of_year)) +
  geom_line(aes(y = RMSE_test), size = 1, colour = "thistle2") +
  geom_line(aes(y = deviance_explained), size = 1, colour = "skyblue1") +
  theme_bw()


