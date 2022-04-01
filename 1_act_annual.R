# to run scripts all scripts (1_act_annaulR - 3_model_eval.R scripts) in terminal:
# 1. change directory: cd ~/Documents/Post\ Doc/Study\ Projects/ACT\ TRAP\ MM/ACT\ HEI\ Supp/act_hei_aim1a
# 2. run bash script: bash run_scripts.bash


# to run individual script in terminal:
# 1. in terminal, change directory: 
#   cd ~/Documents/Post\ Doc/Study\ Projects/ACT\ TRAP\ MM/ACT\ HEI\ Supp/act_hei_aim1a
# 2. enter in terminal: Rscript 1_act_annual.R 

# examples of parallel processing: https://jennybc.github.io/purrr-tutorial/ls12_different-sized-samples.html 

##################################################################################################
# **Purpose** 
#   
#   * We want to better understand how sampling design may impact air pollution predictions models. We use data collected in the ACT TRAP mobile monitoring campaign to conduct sampling simulations and answer questions related to:   
#   - visits per site    
# - stop sampling density 
# - the size of the monitoring area relative to the cohort locations    
# 
# * We noted that prediction models in the California simulations did not always perfom well. Conducting simulations with the ACT TRAP campaign data rather than CA regulatory data should be more revealing since our models perform well


##################################################################################################
# SETUP
##################################################################################################
#tictoc::tic()

# Clear workspace of all objects and unload all extra (non-base) packages
rm(list = ls(all = TRUE))
if (!is.null(sessionInfo()$otherPkgs)) {
  res <- suppressWarnings(
    lapply(paste('package:', names(sessionInfo()$otherPkgs), sep=""),
      detach, character.only=TRUE, unload=TRUE, force=TRUE))
}

pacman::p_load(tidyverse,
               parallel, #mclapply; detectCores()
               future.apply, #future_replicate()
               lubridate # %within%
               )    

set.seed(1)

## for future.apply::future_replicate()
plan(multisession, workers = 6)
 
##################################################################################################
# LOAD DATA
##################################################################################################
#stop data with temporal variables included
stops_all <- readRDS(file.path("..", "..", "2. Prediction Models", "prediction_models_R", "Output", "stops_with_time_vars.rda")) %>%
  #drop duplicate UFP instruments
  filter(!variable %in% c("pmdisc_number", "pnc_screen", "ns_total_conc")) %>%
  #don't need these
  select(-c(contains(c("pollutant", "range")))) %>%
  # add route
  mutate(route = as.numeric(gsub(".*R0", "", runname))) 

# only use non-test sites for simulations
non_test_sites <- readRDS(file.path("Output", "mm_cov_train_set.rda")) %>%
  distinct(location) %>% pull()

stops <- filter(stops_all, location %in% non_test_sites)

##################################################################################################
# COMMON VARIABLES
##################################################################################################
# number of simulations
sim_n <- 30

rush_hours <- c(7:10, 15:18)
business_hours <- c(9:17)


unique_seasons <- unique(stops$season) %>% as.character()
variable_names <- unique(stops$variable)

# number of samples for fewer hours, seasons, reduced balance

fewer_hrs_seasons_n <- 12

##################################################################################################
# ONLY KEEP STOPS W/ ALL POLLUTANTS
##################################################################################################

keep_times <- stops %>%
  distinct(time, variable, value) %>%
  spread(variable, value) %>% #View() #8150 original stops
  mutate(original_stops = n()) %>%
  drop_na() %>%
  mutate(prop_remaining_stops = n()/original_stops) %>% #7345 stops left (90.12%)
  distinct(time) %>% pull()

#range(keep_times) # "2019-03-08 13:51:35 PST" "2020-03-17 23:12:29 PDT"



stops <- filter(stops, time %in% keep_times)

stops_w <- spread(stops, variable, value)

##################################################################################################
# TRUE ANNUAL AVERAGE
##################################################################################################

 
true_annual <- stops_w %>%
  group_by(location, route) %>%
  mutate(visits = n()) %>%
  summarize_at(vars(variable_names, visits), ~mean(.)) %>%
  mutate(
    campaign = 1,
    design = "full",
    version = "all training  data"
  ) %>%
  ungroup()
  
   


# save est set annual averages for validation later
annual_test_set <- stops_all %>%
  filter(!location %in% non_test_sites,
         #same date range as training data
         time >= min(as.Date(keep_times)) & time <= max(as.Date(keep_times))
         ) %>%
  group_by(variable, location) %>%
  summarize(value = mean(value),
            visits = n(),
            campaign = 1,
            design = "test set",
            version = "test set"
  )

saveRDS(annual_test_set, file.path("Output", "annual_test_set.rda"))

##################################################################################################
# SAMPLING DESIGNS
##################################################################################################
# FEWER VISITS

## 1. reduced temporal balance
### A. season-tow2 balanced (not by hour). (12 samples)
### B. season-balanced (not by tow2 or hour). (4 samples)

### --> note: this doesn't consider aggregating samples to 
#              e.g., sample for first 2 wks of the season and then being done


## 2. fewer random visits: seq(5,25,5)

##################################################################################################
## functions 

### sampling functions
s_tow2_sample_fn <- function(df, wkday_visits = 2, wkend_visits = 1) {
  
  no_visits <- ifelse(first(df$tow2) == "Weekday", wkday_visits, wkend_visits)
  df0 <- slice_sample(df, n  = no_visits)
  
  return(df0)
  }

# fn takes one random sample from each list item, according to FUN, and calculates a 'value' average
#list should be in wide format (variable/pollutant names) so that same [temporal] samples are collected across pollutants
                             #list   #sample fn    #descriptor
one_sample_avg <- function(my_list,# = stops_w_list, 
                             my_sampling_fn # =s_tow2_sample_fn  
                             ) {
  result <- mclapply(my_list, mc.cores = 6, FUN = my_sampling_fn) %>%
    #unlist results 
    bind_rows() %>%
    group_by(location, route) %>%
    mutate(visits = n()) %>%
    #calculate annual average
    summarize_at(vars(variable_names, visits), ~mean(.))
    
    return(result)
}

##################################################################################################
#repeat sampling many times

## 1. reduced temporally-balanced sampling 

### A. 12 samples: 4 season x 1 sample/wkend or 2 samples/wkday 
balanced_season_tow2 <- future_replicate(n = sim_n, simplify = F,
                    expr = one_sample_avg(my_list = group_split(stops_w, #variable, 
                                                                location, season, tow2), 
                                            my_sampling_fn = s_tow2_sample_fn)
                    ) %>%
  #unlist
  bind_rows() %>%
  #add simulation number & label
  group_by(location, #variable,  visits
           ) %>%
  mutate(campaign =  row_number(),
         design = "reduced balance",
         version = "balanced by season & TOW2") %>%
  ungroup()
 

# B. 1 sample per season
balanced_season <- future_replicate(n = sim_n, simplify = F,
                           expr = one_sample_avg(my_list = group_split(stops_w, #variable, 
                                                                       location, season), 
                                                 my_sampling_fn = function(x) slice_sample(x, n  = fewer_hrs_seasons_n/4))
                      ) %>%
  #unlist
  bind_rows() %>%
  #add simulation number & label
  group_by(location, #variable,  visits
           ) %>%
  mutate(campaign =  row_number(),
         design = "reduced balance",
         version = "balanced by season")%>%
  ungroup()

##################################################################################################
#2. fewer seasons. keep the number of samples the same. 6 is approximately the # of samples/site/season (i.e., max for season ==1)

## notice that a few sites (e.g., MS0138 in winter) have < 6 samples/season, so we'll sample w/ replacement to keep the numebr of samples the same
# stops_w %>% group_by(location, season) %>% summarize(n = n()) %>% ungroup() %>% summary()

season_n <- c(1:3)

#x = group_split(stops_w, location)[[1]]

season_times <- data.frame()

for (i in seq_along(season_n)) {
  #i=1
  
  temp <- future_replicate(n = sim_n,  
                           simplify = F,
                                    expr = one_sample_avg(my_list = group_split(stops_w, #variable, 
                                                                                location), 
                                                          my_sampling_fn = function(x)  {
                                                            #? keep to make sure locations pick same season combos across simulations ?? or will lapply do this anyways?
                                                            #set.seed(1)
                                                            seasons_to_sample <- sample(c("spring", "summer", "fall", "winter"), size = season_n[i], replace = F)
                                                            
                                                            df <- filter(x, season %in% seasons_to_sample) %>%
                                                              group_by(season) %>%
                                                              # a few sites have < 6 samples/season, so use replacement=True here
                                                              slice_sample(n = fewer_hrs_seasons_n/season_n[i], replace=T)
                                                            }
                                                          )
                           ) %>%
    #unlist
    bind_rows() %>%
    #add simulation number & label
    group_by(location) %>%
    mutate(campaign =  row_number(),
           design = "fewer seasons",
           version = paste0(season_n[i])
           )%>%
    ungroup()
  
  season_times <- rbind(season_times, temp)

}

# OLD: 
# ## 2. shorter campaign duration 
# ### ISSUE: diff # of seasons has diff # visits/site
# season_times0 <- c(list(c("fall", "summer", "winter"), c("spring", "summer", "winter"),
#                         c("summer", "winter"), c("spring", "fall")
#                         ),
#                    unique_seasons)
#   
# season_times <- mclapply(season_times0, mc.cores = 6, FUN = function(x) {
#   stops_w %>% filter(season %in% x) %>%
#     mutate(version = paste0(length(x), ": ", paste0(x, collapse = ", ")),
#            campaign =  1,
#            design = "fewer seasons",
#            )
#   }
#   ) %>%
#   bind_rows() %>%
#   group_by(location, route, campaign, design, version) %>%
#   mutate(visits = n())%>%
#   group_by(location, route, visits, campaign, design, version) %>%
#   summarize_at(vars(variable_names), ~mean(.)) %>%
#   ungroup() %>%
#   select(names(balanced_season))

##################################################################################################
## randomly select fewer samples for the campaign (easy to simulate but less impractical?)
### seq(5, 25, 5)

visit_n <- seq(5, 25, 5)

random_fewer <- list()

for(i in seq_along(visit_n)) {
  #i=1
  temp <- future_replicate(n = sim_n, simplify = F,
                                   expr = one_sample_avg(my_list = group_split(stops_w, location),
                                                         my_sampling_fn = function(x) slice_sample(x, n  = visit_n[i]))
  ) %>%
    #unlist
    bind_rows() %>%
    #add simulation number & label
    group_by(location, visits) %>%
    mutate(campaign =  row_number(),
           design = "fewer visits",
           version = paste0(visit_n[i], "")
           )
  
  random_fewer[[i]] <- temp
  
  }

random_fewer <- bind_rows(random_fewer) %>% ungroup()

##################################################################################################

# only collect same number of samples 


rh_bh <- list(sort(unique(c(rush_hours, business_hours))),
              business_hours, rush_hours
              )

names(rh_bh) <- c("business & rush", "business", "rush")

#x = group_split(stops_w, location)[[1]]

rh_bh_df <- data.frame()

for(i in seq_along(rh_bh)) {
  #i=1
  temp <- future_replicate(n = sim_n,
                           simplify = F,
                           expr = one_sample_avg(my_list =group_split(stops_w, location), 
                           #mc.cores = 4, 
                           my_sampling_fn = function(x) {
                             x %>% filter(tow2 == "Weekday",
                                          hour %in% rh_bh[[i]]) %>%
                               # NOTE: using sampling w/ replacement to ensure 10 samples/site 
                               slice_sample(n=fewer_hrs_seasons_n, replace=T)
                             }
                           )
                           ) %>%
    bind_rows() %>%
    group_by(location) %>%
    mutate(
      campaign = row_number(),
      version = names(rh_bh)[i],
      design = "fewer hours"
    ) %>%
    as.data.frame()
  
  rh_bh_df <- rbind(rh_bh_df, temp)
  
  }







                  
# rh_bh_df <- mclapply(rh_bh, mc.cores = 4, FUN = function(x) {
#                         stops_w %>% filter(tow2 == "Weekday",
#                                            hour %in% x) %>%
#                           mutate(version = ifelse(length(x) == length(rush_hours), "Rush Hours",
#                                                   ifelse(length(x) == length(business_hours), "Business Hours", "Business & Rush Hours")
#                                                   ),
#                                  campaign =  1,
#                                  design = "fewer hours") } ) %>%
#   bind_rows() %>%
#   group_by(location, route, design, version, campaign) %>%
#   mutate(visits = n()) %>%
#   group_by(location, route, visits, design, version, campaign) %>%
#   summarize_at(vars(variable_names), ~mean(.)) %>%
#   ungroup() %>%
#   select(names(balanced_season))



##################################################################################################
# combine TEMPORAL simulation results
temporal_sims <- rbind(
  true_annual,
  balanced_season_tow2,
  balanced_season,
  season_times,
  random_fewer,
  rh_bh_df
  ) %>%
  mutate(spatial_temporal = ifelse(grepl("full", design), "spatial temporal", "temporal"))


##################################################################################################
## SPATIAL simulations: fewer stops (lower density) 
##################################################################################################
### fewer number of stops/sites 

site_n <- c(#10, 
  seq(50, 250, 50)) 

fewer_sites <- future_replicate(n = sim_n, simplify = F, 
                                expr = mclapply(site_n, mc.cores = 6, FUN = function(x) {
                                  slice_sample(true_annual, n  = x) %>%
                                    mutate(version = paste0(x, ""),
                                           design = "fewer sites",
                                           )}) ) %>%
    bind_rows() %>%
    mutate(campaign = rep(1:sim_n, each=sum(site_n)))
      
##################################################################################################
## smaller area/fewer routes
### many sources: 1 (urban center, many participants), 6 (airport, Tacoma/suburban/industrial/port)   

### --> using a loop in case we add more route combinations

routes_list <- list(c(1:2), c(1:4), c(1:6)) #list(c(1:3), c(1:6))

fewer_routes <- mclapply(routes_list, mc.cores = 6, FUN = function(x) {
  filter(true_annual, route %in% x) %>%
    mutate(version = paste(range(x), collapse = "-"),
           design = "fewer routes",
           campaign = 1)} ) %>%
  bind_rows()  

##################################################################################################
# combine spatial simulations
spatial_sims <- rbind(fewer_sites, fewer_routes) %>%
  mutate(spatial_temporal = "spatial")

##################################################################################################
# combine spatial and temporal simulations

annual_training_set <- rbind(temporal_sims, spatial_sims)

##################################################################################################
# SAVE DATA
##################################################################################################
saveRDS(annual_training_set, file.path("Output", "annual_training_set.rda"))

#tictoc::toc()

print("done with 1_act_annual.R")


