
# Plot simulation results

# Script show how plots in "Contouring Groundwater Information: a Search for Truth"
#   are coded. All figures generated in this script are examples of how the figures
#   are generated but do not directly recreate exact Figures in the manuscript

# libraries --------------------------------------------------------------------
library(tidyverse)
library(viridis)

rm(list=ls())

# load data --------------------------------------------------------------------
plume_type <- "complex"

result_dfs <- list.files(
  path = "data",
  pattern = paste0("sim_", plume_type,"_df.*\\.RDS$"),
  full.names = TRUE,
  recursive = TRUE
)

load_and_get_data <- function(filepath) {
  tryCatch({
    # create a temporary environment
    temp_env <- new.env()
    
    # load the .RDS file into the temporary environment
    load(filepath, envir = temp_env)
    
    # get the names of the objects that were loaded
    loaded_object_names <- ls(envir = temp_env)
    
    if (length(loaded_object_names) == 0) {
      warning(paste("No objects loaded from file:", filepath))
      return(NULL)
    }
    
    # extract the objects from the temporary environment
    loaded_objects <- mget(loaded_object_names, envir = temp_env)
    
    # if there's only one object, return it directly; otherwise, return the list.
    if (length(loaded_object_names) == 1) {
      return(loaded_objects[[1]])
    } else {
      return(loaded_objects)
    }
    
  }, error = function(e) {
    message(paste("Error loading file:", filepath, "\n", e))
    return(NULL) # Return NULL on error
  },
  warning = function(w){
    message(paste("Warning loading file", filepath, "\n",w))
    return(readRDS(filepath))
  })
}

# Load all .RDS files and store in a list
results_list <- lapply(result_dfs, load_and_get_data)
mspe_all <- do.call(rbind, results_list)

# Figure 9: TPM grid only boxplots ---------------------------------------------
mspe_all %>%
  filter(full_grid == FALSE) %>%
  rename(Method = "pred") %>%
  mutate(Method = case_when(
    Method == "rCCR" ~ "rCCR",
    Method == "rsfi" ~ "RFSI",
    Method == "tps_time" ~ "TPS-T",
    Method == "SB" ~ "BAYES", 
    .default = toupper(Method)
  )) %>%
  ggplot() +
  geom_boxplot(aes(y = mspe, x = Method, fill = Method),
               outliers = FALSE) +
  # scale_fill_viridis(discrete = T) +
  scale_fill_manual(values = palette.colors(palette = "Okabe-Ito")) +
  facet_wrap(vars(n_wells),
             scales = "free") +
  scale_y_continuous(trans='log10') +
  labs(y = "MSPE") +
  theme_bw()

# Table of percent of simulations that outperform the others -------------------

perc_best <- mspe_all |> 
  filter(time %in% c("ts10", "ts20", "ts30", "ts40")) |> # keep only 4 times points
  group_by(full_grid, network_id, n_wells, time) |> 
  slice_min(mspe) |> # take minimum mspe for each group
  ungroup() |> 
  group_by(full_grid, n_wells, time) |> 
  summarise(
    tpm_best = mean(pred == "tpm"),
    tps_best = mean(pred == "tps"),
    tpst_best = mean(pred == "tps_time"),
    sb_best = mean(pred == "SB"),
    rccr_best = mean(pred == "rCCR"),
    inla_best = mean(pred == "inla"),
    rfsi_best = mean(pred == "rsfi")
  ) |> 
  mutate(
    across(ends_with("best"), function(x) round(x * 100, 1)),
    tpm_best = case_when(full_grid ~ "", !full_grid ~ as.character(tpm_best))
  ) |> 
  arrange(full_grid, -n_wells) |> 
  ungroup() |> 
  mutate(
    full_grid = case_when(
      full_grid == FALSE & time == "ts10" & n_wells == 100 ~ "TPM Region",
      full_grid == TRUE & time == "ts10" & n_wells == 100 ~ "Full Grid",
      TRUE ~ ""
    )
  )

perc_best <- perc_best |> 
  mutate(n_wells = case_when(time != "ts10" ~ "", TRUE ~ as.character(n_wells))) |> 
  mutate(time = str_extract(time, "\\d+")) |> 
  rename(
    "Prediction Region" = full_grid,
    "No. of Wells" = n_wells,
    "Time" = time,
    "TPM" = tpm_best,
    "TPS" = tps_best,
    "TPS-T" = tpst_best,
    "BAYES" = sb_best,
    "rCCR-S" = rccr_best,
    "INLA" = inla_best,
    "RFSI" = rfsi_best
  )  

view(perc_best)


# Figure 3: contour example from each method -----------------------------------
# load "true" data
sim_file_path <- "data"
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

z_lim <- range(plume_long$conc)

# load result data -- manuscript has two sets of results, only one is plotted
#   here as an example
load("data/")
str(pred_df_out)

# subset to 10, 20, 30, 40 ts
pred_subset <- pred_df_out %>%
  filter(time %in% c("ts10", "ts20", "ts30", "ts40"))
saveRDS(pred_subset, file = "complex_pred_plume_network_25_100_nwells_subset.RDS")


# select number of wells -- can be 15, 30, 50 75, or 100
n_wells <- 100

# plot
pred_df_out %>%
  filter(time %in% c("ts10","ts20","ts30", "ts40")) %>%
  rename(Method = "pred", Concentration = "conc") %>%
  mutate(Method = case_when(
    Method == "rCCR" ~ "rCCR",
    Method == "rsfi" ~ "RFSI",
    Method == "tps_time" ~ "TPS-T",
    Method == "SB" ~ "BAYES",
    .default = toupper(Method)
  )) %>%
  ggplot() +
  geom_tile(data = plume_long %>% 
              filter(time %in% c("ts10","ts20","ts30", "ts40")) %>%
              rename(Concentration = "conc"), 
            aes( x = x, y = y, fill = Concentration )) +
  geom_contour(data = plume_long %>% 
                 filter(time %in% c("ts10","ts20","ts30", "ts40")) %>%
                 rename(Concentration = "conc"),
               aes(x = x, y = y, z = Concentration),
               breaks = 20, linewidth = 0.8, col = "white") +
  scale_fill_gradient(low = "lightyellow", high = "darkblue", limits = z_lim) +
  # scale_fill_gradientn(colours = fields::tim.colors(10)) +
  geom_contour(aes(x = x, y = y, z = Concentration, col = Method
                   # linetype = Method
                   ),
               breaks = 20, linewidth = 0.8) +
  # scale_color_viridis(discrete = T) +
  scale_color_manual(values = palette.colors(palette = "Okabe-Ito")) +
  facet_wrap(vars(time)) +
  coord_equal() +
  lims(y = c(250, 800)) +
  labs(title =paste0("Simulated Well Contours with ", n_wells," wells")) +
  theme_bw() +
  theme(legend.position = "bottom")



# Figures 4-8: Contour spatial boxplots example --------------------------------
# example data is available for 100 wells only

# load example result data directly from data folder
df <- load_and_get_data("gw_contours/data/predictions_complex_plume.RDS")

# get average, median, 5th, 95th contours by pixel -- these lines make take a minute to run
df_fbplot <- df %>% 
  group_by(x, y, time, pred) %>%
  summarise(
    # conc_avg = mean(conc, na.rm = TRUE),
    conc_med = median(conc, na.rm = TRUE),
    conc_05 = quantile(conc, probs = c(0.05), na.rm = TRUE),
    conc_95 = quantile(conc, probs = 0.95, na.rm = TRUE),
    conc_se = sd(conc)
  ) %>%
  ungroup()

# get true data from df and store as its on data frame
df_true <- df %>% 
  dplyr::select(x, y, obs_conc, time) %>% 
  distinct()

# range of true data will be used for plot range
z_lim <- range(df_true$obs_conc)

# set the color palette
cols <- palette.colors(palette = "Okabe-Ito")

# plot -- note: true contour is on the plot but not included in legend
df_fbplot %>%
  mutate(pred = case_when(
    pred == "rCCR" ~ "rCCR",
    pred == "rsfi" ~ "RFSI",
    pred == "tps_time" ~ "TPS-T",
    pred == "SB" ~ "BAYES",
    .default = toupper(pred)
  )) %>%
  ggplot() +
  geom_tile(data = df_true,
            aes(x = x, y = y, fill = obs_conc)) +
  geom_contour(data = df_true,
               aes(x = x, y = y, z = obs_conc),
               col = "white",
               breaks = 20) +
  geom_contour(aes(x = x, y = y, z = conc_med, color = "Median"), 
               breaks = 20) +
  geom_contour(aes(x = x, y = y, z = conc_05, color = "5th Percentile"), 
               breaks = 20) +
  geom_contour(aes(x = x, y = y, z = conc_95, color = "95th Percentile"), 
               breaks = 20) +
  facet_grid(rows = vars(pred),
             cols = vars(time)) +
  scale_color_manual(name = "Contours", # Set the legend title
                     values = c("Median" = cols[1],
                                "5th Percentile" = cols[2],
                                "95th Percentile" = cols[7]
                     )) +
  scale_fill_gradient(low = "lightyellow", high = "darkblue", limits = z_lim,
                      name = "Concentration") +
  ylim(c(0, 1000)) +
  labs(x = "Column", y = "Row") +
  coord_equal() +
  theme_bw()

  

  







