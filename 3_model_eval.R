#script purpose: evaluate UK-PLS model performances

##################################################################################################
# to run in terminal:
# 1. in terminal, change directory: 
#   cd ~/Documents/Post\ Doc/Study\ Projects/ACT\ TRAP\ MM/ACT\ HEI\ Supp/act_hei_aim1a
# 2. enter in terminal: Rscript 3_model_eval.R 

##################################################################################################
# SETUP
##################################################################################################
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
# Lambert Conic projection (meters)
lambert_proj <- "+proj=lcc +lat_1=33 +lat_2=45 +lat_0=39 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"

predictions <- readRDS(file.path("Output", "predictions.rda"))

















##################################################################################################
# CV STATS FUNCTION
##################################################################################################
# # Define function to create a bubble plot for kriging residuals
# krige.cv.bubble <- function(cv.out, plot_title){
#   ggplot(data=cv.out) +
#     geom_sf(aes(size=abs(residual), color=factor(residual>0)), alpha=1/2) +
#     scale_color_discrete(name='residual > 0', direction=-1) +
#     scale_size_continuous(name='|residual|') +
#     ggtitle(plot_title) +
#     theme_bw()
# }

### --> UPDAT FN 

# Define function for calculating the MSE, RMSE, MSE-based R2 from krige.cv output
krige.cv.stats <- function(krige.cv.output, description){
  d <- krige.cv.output
  # mean of observations
  mean_observed <- mean(d$observed)
  # MSE of predictions
  MSE_pred <- mean((d$observed - d$var1.pred)^2)
  # MSE of observations (for R2 denominator)
  MSE_obs <- mean((d$observed - mean_observed)^2)
  # # print the results not rounded
  # cat(paste("RMSE:  ", sqrt(MSE_pred)),'\n')
  # cat(paste("MSE-based R2:  ", max(1 - MSE_pred/MSE_obs, 0)),'\n')
  
  # make dataframe of performance results
  tibble(Description = description, 
         RMSE = round(sqrt(MSE_pred), 4), 
         MSE_based_R2 = round(max(1 - MSE_pred/MSE_obs, 0), 4) )
}



##################################################################################################

# calculate performance statistics across CVs
krige.cv.stats(dt.CV5uk, "UK on sqrt(dist): 5-fold CV") 
##################################################################################################












##################################################################################################
# CALCULATE GEOGRAPHIC DISTANCES
##################################################################################################

# o	Additionally, for the spatial sims
# 	Fewer random sites sims: RMSE vs Monitoring density (monitors in monitoring area/km2); avg distance between training monitors; avg distance between training and test (“cohort”) monitors 



