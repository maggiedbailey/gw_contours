# apply each method to simulated data set

# clear environment
rm(list=ls())

options(dplyr.summarise.inform = FALSE)

# load libraries
library(MASS)
library(gstat)
library(whitening)
library(sp)
library(fields)
library(interp)
library(tictoc)
library(meteo)
library(tidyr)
library(tidyverse); theme_set(theme_minimal())

# setwd("github_repo/")

# source functions -------------------------------------------------------------
source("model_functions.R")
source("inla_functions.R") # model functions depends on funcs in this script

# read in simulation data ------------------------------------------------------
sim_file_path <- "data/"
plume_type <- "complex" # complex | mid | simple
tmp_list <- lapply(
  grep(plume_type,
       list.files( sim_file_path, full.names = TRUE, pattern = "*.rds"),
       value = TRUE),
  readRDS
  )
plume_long <- do.call(rbind, tmp_list)

rm(tmp_list)

# add time_dbl
plume_long <- plume_long %>%
  mutate(
    time_dbl = as.numeric(ts), # add time_dbl as number for time step
    time = paste0("ts", ts) # add time as character string for time step
    ) %>%
  dplyr::select(-ts)

# one row per well
plume_wide <- plume_long |>
  pivot_wider(
    names_from = "time",
    values_from = "conc",
    id_cols = c(x,y)
  )

# check in the correct orientation
plume_long |>
  ggplot() +
  geom_raster(aes(x = x, y = y, fill = conc)) +
  scale_fill_viridis_c() +
  facet_wrap(vars(time_dbl)) +
  theme(panel.grid = element_blank())

# all grid locations - data frame with two columns
plume_xy <- plume_wide |> dplyr::select(x, y) 


# run simulation ---------------------------------------------------------------

# number of wells in the domain to use as observations
well_interest <- 15 # 15, 30, 50 75, 100 but set 15 here for example

k <- 1 # number of well networks to simulate
d <- 1:k
x <- seq_along(d)
max_n <- 8 # 8 networks in each rds
iter_i <- split(d, ceiling(x/max_n))
iter_n_names <- 1:length(iter_i) # number of iteration groups with 8 networks

for (j in iter_n_names) {
  
  cat(" --- iteration group ", j, " ---", fill = TRUE)

    tic("group", as.character(j))
    
    # creating groups of 8 for jth group
    vec <- (j * 8 - 7):(j * 8) 
    
    # where to store results for contouring r&d
    file_path <- file.path("data")

    for (k in vec) { 
      
      cat("k: ", k, fill = TRUE)
      sim_df <- simulate_well_network(
        k,
        well_n_int = well_interest,
        plume_data_long = plume_long,
        plume_type = plume_type,
        plume_data_wide = plume_wide
      )
      obj_name <- paste0("sim_", plume_type, "_df_", k,"_itergroup_", j)
      assign(obj_name, sim_df)
      save(list = obj_name, file = paste0(file_path, "/", obj_name, ".RDS"))

    }
    
    toc()
}









  