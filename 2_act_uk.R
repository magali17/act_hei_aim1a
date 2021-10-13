# Script Purpose:
# * built UK models & predict at out-of-sample locations

##################################################################################################
# to run in terminal:
# 1. in terminal, change directory: 
#   cd ~/Documents/Post\ Doc/Study\ Projects/ACT\ TRAP\ MM/ACT\ HEI\ Supp/act_hei_aim1a
# 2. enter in terminal: Rscript 2_act_uk.R 

##################################################################################################
# NOTES
##################################################################################################

### --> save & drop code for clustering if don't use


### --> ??
# **Note on `rgdal`**:  If you are running your lab locally on a Mac, the `rgdal`
# package doesn't load properly on newer Macs running Catalina.  You need to
# install GDAL from the terminal.  See this
# [link](https://medium.com/@egiron/how-to-install-gdal-and-qgis-on-macos-catalina-ca690dca4f91)
# for how to do this.

 
  


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
               #future.apply, #future_replicate()
               pls, 
               gstat, #variogram()
               sf#, #for spatial data; st_
               #sp #has muse datasedt       # ? nEED?
)    

set.seed(1)

# ## for future.apply::future_replicate()
# plan(multisession, workers = 6)

##################################################################################################
# LOAD DATA
##################################################################################################
# mapping variables
project_crs <- 4326  #lat/long
m_crs <- 32148

# save coordinate systems as variables
# WGS84 latitude-longitude
#latlong_proj <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
# Lambert Conic projection (meters)
lambert_proj <- "+proj=lcc +lat_1=33 +lat_2=45 +lat_0=39 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"

##################################################################################################
# mm annual estimates
annual <- readRDS(file.path("Output", "annual_training_set.rda")) %>%
  #add covariates
  left_join(readRDS(file.path("Output", "mm_cov_train_set.rda"))) %>%
  #convert to sf
  st_as_sf(coords = c('longitude', 'latitude'), crs=project_crs, remove = F) %>%
  st_transform(m_crs) %>%
  #log-transform UFPs since these are right skewed. all other variables look fine.
  #mutate(pnc_noscreen = log(pnc_noscreen)) %>%
  gather("variable", "value", ma200_ir_bc1:pnc_noscreen) %>%
  relocate(variable, value, .before = longitude)

## same for test set
annual_test_set <- readRDS(file.path("Output", "annual_test_set.rda")) %>%
  #add covariates
  left_join(readRDS(file.path("Output", "mm_cov_test_set.rda"))) %>%
  #convert to sf
  st_as_sf(coords = c('longitude', 'latitude'), crs=project_crs, remove = F) %>%
  st_transform(m_crs) #%>%
  #log-transform UFPs since these are right skewed. all other variables look fine.
  #mutate(value = ifelse(variable == "pnc_noscreen", log(value), value))
 
##################################################################################################
# COMMON VARIABLES
##################################################################################################
  
# for modeling 
cov_names <- st_drop_geometry(annual) %>% ungroup() %>%
  select(log_m_to_a1:last_col()) %>% names() # 224 covarites

pls_comp_n <- 2

#k-folds for CV
k=10

var_names <- unique(annual$variable)

##################################################################################################
# SETUP
##################################################################################################

# # everything look fairly normally distributed after UFPs are log-transformed 
# annual %>%
#   filter(design == "full training campaign") %>%
#   ggplot(aes(x=value)) +
#   facet_wrap(~variable, scales="free") +
#   geom_histogram() +
#   labs(title = "Distribution of site annual average concentrations")

 

##################################################################################################
# CREATE VALIDATION CLUSTERS FOR EACH VARIABLE-DESIGN-VERSION
##################################################################################################

# spatial_cluster <- function(dt, k=k) {
#   
#   # create spatial clusters with the available data for validation 
#   long_lat <- c("longitude", "latitude")
#   
#   site_lat_log <- dt %>% st_drop_geometry() %>% ungroup() %>%
#     #variables used for clustering
#     distinct(location, longitude, latitude,
#              spatial_temporal, design, version, campaign
#              ) %>%
#     # may not need to scale, but groups look slightly better.
#     mutate_at(long_lat, ~scale(.))
#   
#   #if sites arranged in the same order & simulations have the same sites (i.e., temporal simulations), 
#   # sites will end up in same folds across simulations/designs
#   set.seed(2)
#   
#   # using Lloyd algorithm, which allows for k=nrow(df), otherwise 10 visits has issues  
#   clusters <- kmeans(select(site_lat_log, long_lat), centers = k, algorithm = "Lloyd", 
#                      #increase from default of 10. otherwise does not always converge.
#                      iter.max = 30)
#   
#   site_lat_log <- site_lat_log %>%
#     mutate(cluster = clusters$cluster) %>%
#     select(-long_lat)
# 
#   return(site_lat_log)
#   
#   # print("number of sites in each cluster:")
#   # summary(clusters$size)
#   # 
#   # dt %>%
#   #   mutate(cluster = factor(cluster)) %>%
#   #   ggplot(aes(col=cluster)) +
#   #   geom_sf() +
#   #   labs(title = "Example of clustered validation groups for 'full campaign'")
#   # 
#   # ggsave(file.path("Images", "clustered_cv.png"), width = 6, height = 9)
# 
# }


############################################################################################################

# ### --> ? only do this for specific spatial simulation? 
# 
# cluster_df <- lapply(group_split(annual, spatial_temporal, design, version, campaign),
#                             function(x) spatial_cluster(x, k=k)) %>%
#   bind_rows()
#   
# #join clusters to annual 
# annual <- suppressMessages(left_join(annual, cluster_df)) %>%
#   select(cluster, everything())
 
##################################################################################################
# CREATE RANDOM VALIDATION FOLDS FOR EACH VARIABLE-DESIGN-VERSION
##################################################################################################

random_fold <- function(df, k.=k) {
  #make sure temporal sims receive same fold designation
  set.seed(2)
  
  result <- df %>% st_drop_geometry() %>%
    distinct(location,
             spatial_temporal, design, version, campaign
             ) %>%
    mutate(random_fold = sample(1:k.,size = nrow(.), replace = T, ))
  
  return(result)
  }

############################################################################################################

### --> ? only do this for specific spatial simulation? 

random_fold_df <- lapply(group_split(annual, spatial_temporal, design, version, campaign),
                     function(x) random_fold(x, k.=k)) %>%
  bind_rows()

#join clusters to annual 
annual <- suppressMessages(left_join(annual, random_fold_df)) %>%
  select(#cluster, 
         random_fold, everything())


############################################################################################################
# UK-PLS FUNCTION 
############################################################################################################
# fn returns UK-PLS predictions. inputs are two spatial objects (simple features). fn automatically transforms these to a lambert projection

 
# modeling_data = x[[1]]
# new_data = filter(annual_test_set, variable == var_names[i])

uk_pls <- function(modeling_data, # data for fitting pls-uk models
                   new_data, #prediction locations
                   cov_names. = cov_names,  #covariates to be used in modeling
                   pls_comp_n. = pls_comp_n #pls components to use
                   ) {
  
  #lambert projection for UK model
  lambert_proj <- "+proj=lcc +lat_1=33 +lat_2=45 +lat_0=39 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
   ############################################################################################################
  # fit PLS model to estimate fewer components from geocovariates
  
  print(paste(first(modeling_data$version), first(modeling_data$variable), "campaign ", first(modeling_data$campaign)#, "test set ", c(1:10)[!c(1:10) %in% unique(modeling_data$cluster)]
              ))
  
  pls_model <- plsr(as.formula(paste('value ~', paste(cov_names., collapse = "+"))),
                    data = modeling_data, ncomp = pls_comp_n., scale=T, center=T)
  
  #extract compoent scores for UK
  modeling_data_scores <- predict(pls_model, type = "scores") %>% data.frame() %>% #head()
    # add location & value info
    cbind(data.frame(select(modeling_data, -cov_names.))) %>%
    #convert back to sf. geom is dropped otherwise
    st_as_sf()
  
  new_data_scores <- predict(pls_model, type = "scores", newdata = new_data) %>% data.frame() %>%  
    # add location & value info
    cbind(data.frame(select(new_data, -cov_names.))) %>%
    #convert back to sf. geom is dropped otherwise
    st_as_sf()
  
  ############################################################################################################
  # fit UK models & predict at new locations
  
  # UK formula 
  uk_formula <- as.formula(paste("value ~", paste0("Comp.", 1:pls_comp_n., collapse = "+")))

  # estimate the variogram model: fit a variogram model, offering to the function several different model options (exponential, spherical, and Matern):
  # using lambert coordinates b/c vertical/horizontal units represent the same jump
  # the default distance in gstat is 1/3 of the maximum distance (use cutoff option to change this)
  
  v.uk <- variogram(uk_formula, st_transform(modeling_data_scores, lambert_proj) )
  m.uk <- fit.variogram(v.uk, vgm(c("Exp", "Sph", "Mat")) )
  #make sure Exp/Sph range estimate is at least 0 when little/no correlation in the data 
  m.uk$range[2] <- max(m.uk$range[2], 1)
  
  #plot(v.uk, m.uk)
  
  # fit UK to the modeling data and predict at the new data locations
  uk_model <- krige(formula = uk_formula, st_transform(modeling_data_scores, lambert_proj), 
                    newdata = st_transform(new_data_scores, lambert_proj), 
                    model = m.uk)
   
  #save predictions
  result <- select(new_data, -cov_names.) %>%
    mutate(prediction = uk_model$var1.pred)
  
  return(result)
  
}
 
##################################################################################################
# CV function
##################################################################################################
# function returns cross-valited predictions for a given dataset

# x = filter(annual, grepl("full", design), grepl("ma200", variable))
# fold_name = "random_fold"

do_cv <- function (x, fold_name) {
  
  k = length(unique(x[[fold_name]]))
  
  df <- data.frame()
  
  for(f in 1:k) {
    #f=1
    modeling_data0 = filter(x, !!as.symbol(fold_name) != f)
    new_data0 = filter(x, !!as.symbol(fold_name) == f)
    temp <- uk_pls(modeling_data = modeling_data0, new_data = new_data0) %>% st_drop_geometry() 
    df <- rbind(df, temp)
  }
  
  return(df)
}

##################################################################################################
# PREDICT
##################################################################################################

common_vars <- c(#"cluster", 
                 "location", "route", "visits", "campaign", "design", "version", "spatial_temporal", "variable", #"value",
                 "prediction"
                 )

# # 1. clustered CV predictions
# ## each df in the list is a simulation; running 10-FCV within each simulation to get CV predictions
# cluster_predictions <- mclapply(group_split(annual, spatial_temporal, design, version, campaign, variable), 
#                       mc.cores = 5, FUN = do_cv) %>%
#   bind_rows() %>%
#   rename(cluster_prediction = prediction)


# 0. random folds
## each df in the list is a simulation; running 10-FCV within each simulation to get CV predictions
random_fold_predictions <- mclapply(group_split(annual, spatial_temporal, design, version, campaign, variable),
                      mc.cores = 5, FUN = do_cv, fold_name = "random_fold") %>%
  bind_rows() %>% 
  select(common_vars) %>%
  mutate(out_of_sample = "Random 10-Fold")  




# 1. test set predictions
 
# WARNNINGS: In fit.variogram(object, x, fit.sills = fit.sills, fit.ranges = fit.ranges,  ... : singular model in variogram fit
# see Note in help(fit.variogram). This has to do with flat variograms w/ little spatial correlation. Try plotting the variograms.

test_set_predictions <- data.frame()

#1 pollutant at a time to make sure things are arranged correctly with the test set
for(i in seq_along(var_names)) {
  
   df <- mclapply(group_split(filter(annual, variable == var_names[i]), spatial_temporal, design, version, campaign), 
                  mc.cores = 5,
                function(x) {
                  df = uk_pls(modeling_data = x, new_data = filter(annual_test_set, variable == var_names[i])) %>%
                    #fn has binding issues later if don't drop geom?
                    st_drop_geometry() %>%
                    #add info to new dataset about the prediction model
                    mutate(
                      spatial_temporal = first(x$spatial_temporal), 
                      design = first(x$design), 
                      version  = first(x$version), 
                      campaign = first(x$campaign)
                    )  
                  }) %>%
     bind_rows()
   
   test_set_predictions <- bind_rows(test_set_predictions, df)
   
   }

test_set_predictions <- test_set_predictions %>% 
  select(common_vars) %>%
  mutate(out_of_sample = "Test")


#2. predict at routes 7-9

route_predictions <- mclapply(group_split(filter(annual, design == "fewer routes"), version, variable), 
                              mc.cores = 5, function(x) {

  non_training_sites <- annual %>%
    filter(variable == first(x$variable),
           route %in% c(7:9),
           #only need one row per site
           grepl("full", design)
           )
  
  uk_pls(modeling_data = x, new_data = non_training_sites) %>%
    st_drop_geometry() %>%
    #add info to new dataset about the prediction model
    mutate(
      spatial_temporal = first(x$spatial_temporal),
      design = first(x$design),
      version  = first(x$version)
      )
  }) %>%
  bind_rows()  

route_predictions <- route_predictions %>%
  select(common_vars) %>%
  mutate(out_of_sample = "Routes 7-9")

 
# # 3. predictions at other, non-training sites
#  
# annual_spatial_list <- annual %>%
#   filter(spatial_temporal == "spatial") %>%
#   group_split(spatial_temporal, design, version, campaign, variable)
#  
# other_sites_predictions <- mclapply(annual_spatial_list, mc.cores = 5, function(x) {
#   training_sites <- unique(x$location)
#   
#   non_training_sites <- annual %>% 
#     filter(!location %in% training_sites,
#            #only need one row per site
#            grepl("full", design), 
#            variable == first(variable))
#   
#   df = uk_pls(modeling_data = x, new_data = non_training_sites) %>%
#     st_drop_geometry() %>%
#     #add info to new dataset about the prediction model
#     mutate(
#       variable = first(x$variable),
#       spatial_temporal = first(x$spatial_temporal), 
#       design = first(x$design), 
#       version  = first(x$version), 
#       campaign = first(x$campaign)
#     ) 
#   }) %>%
#   bind_rows() %>%
#   rename(other_sites_prediction = prediction)



##################################################################################################
# COMBINE PREDICTIONS
##################################################################################################
 
# add campaign-specific and gold-standard estimates to predictions 
# in UK, can't predict at same places, so can't compare estimates & predictions w/o doing CV

#GS estimates for 278 sites
annual_gs_estimates <- annual %>% st_drop_geometry() %>%
  filter(grepl("full", design)) %>%
  distinct(location, variable, value) %>%
  rename(gs_estimate = value)

#GS estimates for 31 sites
test_gs_estimates <- annual_test_set %>% st_drop_geometry() %>%
  distinct(location, variable, value) %>%
  rename(gs_estimate = value)

#GS estimates for all 309 sites
gs_estimates <- bind_rows(annual_gs_estimates, test_gs_estimates)

# estimates from specific campaign simultaions (n=278 sites)
campaign_estimates <- annual %>% st_drop_geometry() %>%
  distinct(location, design, version, campaign, variable, value) %>%
  rename(campaign_estimate = value)


# combine everything
predictions <- rbind(random_fold_predictions, test_set_predictions) %>%
  rbind(route_predictions) %>%
  #left join b/c locations w/ predictions may be fewer than the 309 sites if dno't do 10FCV
  left_join(gs_estimates) %>%
  left_join(campaign_estimates)
  

##################################################################################################
# SAVE DATA
##################################################################################################
saveRDS(predictions, file.path("Output", "predictions.rda"))

#tictoc::toc()
