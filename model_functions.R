# Scripts codes functions for predicting each plume type using the described methods

# method functions -------------------------------------------------------------

# three point method
tpm_pred_function <- function(well_xy, well_long){
  x <- well_xy$x; y <- well_xy$y
  xo <- plume_xy$x[plume_xy$x >= min(x) & plume_xy$x <= max(x)] |> unique()
  yo <- plume_xy$y[plume_xy$y >= min(y) & plume_xy$y <= max(y)] |> unique()
  tpm_pred_df <- plume_xy
  
  tm_interest <- paste0("ts", 1:40)
  
  for (tm in tm_interest){ # tm = tm_interest[1]
    well_time <- filter(well_long, time == tm)
    
    il <-interp(x = well_time$x, y = well_time$y, z = well_time$conc,
                xo = xo, yo = yo)
    idf <- interp2xyz(il, data.frame = TRUE) |> 
      as_tibble() |> 
      rename(!!sym(tm) := z)
    
    tpm_pred_df <- tpm_pred_df |> 
      left_join(idf, by = c("x", "y"))
  }
  
  list(
    "df" = tpm_pred_df,
    "long" = tpm_pred_df |> 
      pivot_longer(starts_with("ts"), names_to = "time", values_to = "conc") |> 
      mutate(pred = "tpm")
  )
}

# tps (no time)
tps_pred_function <- function(well_xy, well_long){
  tps_pred_df <- plume_xy
  
  tm_interest <- paste0("ts", 1:40)
  #tm_interest <- "ts40"
  
  for (tm in tm_interest){
    well_time <- filter(well_long, time==tm)
    mod_tps <- Tps(x = select(well_time, x, y), 
                   Y = well_time$conc)
    pred <- predict(mod_tps, x = plume_xy)[,1]
    tps_pred_df <- tps_pred_df |> 
      mutate(!!sym(tm) := pred)
  }
  
  list(
    "df" = tps_pred_df,
    "long" = tps_pred_df |> 
      pivot_longer(starts_with("ts"), names_to = "time", values_to = "conc") |> 
      mutate(pred = "tps")
  )
}

# tps with time
tps_time_pred_function <- function(well_xy, well_long){
  
  tic()
  mod_tps_time <- Tps(x = well_long[,c("x", "y", "time_dbl")],
                      Y = well_long["conc"])
  toc()

  # create prediction grid
  times <- sort(unique(well_long$time_dbl))
  pred_grid <- expand.grid(
    x = unique(plume_xy$x),
    y = unique(plume_xy$y),
    times = times
  )

  # predict for all time points and years -- this takes a while
  tic("predict")
  out <- predict(mod_tps_time, x = pred_grid)
  toc()

  tps_time_pred_df <- cbind(pred_grid, out)
  tps_time_pred_df <- tps_time_pred_df |>
    mutate(time = paste0("ts", times)) |>
    select(x, y, conc, time)
  
  tps_time_pred_wide <- tps_time_pred_df |>
    pivot_wider(
      id_cols = c(x,y),
      names_from = time,
      values_from = conc
    )
  
  
  list(
    "df" = tps_time_pred_wide,
    "long" = tps_time_pred_df |>
      mutate(pred = "tps_time")
  )
}

# regression model with spatial smoothing
rccr_pred_function <- function(well_xy, well_long){
  cca_prior <- plume_xy
  
  for (tm in paste0("ts", 1:5)){
    well_time <- filter(well_long, time==tm)
    mod_tps <- Tps(x = select(well_time, x, y), 
                   Y = well_time$conc)
    pred <- predict(mod_tps, x = plume_xy)[,1]
    cca_prior <- cca_prior |> 
      mutate(!!sym(tm) := pred)
  }
  
  cca_prior_long <- cca_prior |> 
    pivot_longer(starts_with("ts"), names_to = "time", values_to = "conc") |> 
    mutate(time_dbl = str_extract(time, "\\d+") |> as.numeric()) |> 
    mutate(conc = ifelse(conc < 0, 0, conc))
  
  for (i in 6:40){
    #print(i)
    #make matrices
    x_train <- well_long |> 
      select(x, y, conc, time_dbl) |> 
      filter(time_dbl < i) |> 
      arrange(x, y) |> 
      pivot_wider(names_from = c("x", "y"), values_from = "conc") |> 
      select(-time_dbl) |> 
      as.matrix()
    
    x_test <- well_long |> 
      select(x, y, conc, time_dbl) |> 
      filter(time_dbl == i) |> 
      arrange(x, y) |> 
      pivot_wider(names_from = c("x", "y"), values_from = "conc") |> 
      select(-time_dbl) |> 
      as.matrix()
    
    all_train <- cca_prior_long |> 
      select(-time) |> 
      filter(time_dbl < i) |> 
      arrange(x, y) |> 
      pivot_wider(names_from = c("x", "y"), values_from = "conc")  |> 
      select(-time_dbl) |> 
      as.matrix()
    
    all_train <- all_train + matrix(
      rnorm(length(all_train), 0, sd = 0.000000000001), 
      nrow = nrow(all_train))
    
    #find cca
    cca_mod_all <- scca(x_train, all_train, scale = FALSE, verbose = FALSE)
    cca_mod_all$lambda <- zapsmall(cca_mod_all$lambda) # remove very small values
    cca_pred <- x_test %*% cca_mod_all$PhiX %*% cca_mod_all$K 
    
    cca_mod_pred <- tibble(
      cca_pred_c = cca_pred[1,] |> unname()
    ) |> 
      mutate(
        x = colnames(cca_pred) |> str_extract("[^_]+") |> as.numeric(),
        y = colnames(cca_pred) |> str_extract("(?<=_)\\d+") |> as.numeric(),
      ) |> 
      select(x, y, cca_pred = cca_pred_c)
    
    cca_mod_pred_df <- cca_mod_pred|> 
      right_join(filter(well_long, time_dbl == i), by = c("x", "y")) |> 
      mutate(resid = conc - cca_pred)
    
    #find spatial residual
    resid_tps <- Tps(x = select(cca_mod_pred_df, x, y), 
                     Y = select(cca_mod_pred_df, resid))
    p_r <- plume_xy |> 
      mutate(pred_resid = predict(resid_tps, x = plume_xy)[,1] |> unname()) 
    
    
    #pred df
    pred_df_cca <- cca_mod_pred |> 
      full_join(p_r, by = c("x", "y")) |> 
      mutate(
        conc = cca_pred + pred_resid, 
        conc = ifelse(conc < 0, 0, conc),
        time_dbl = i, 
        time = paste0("ts", time_dbl)
      ) |> 
      select(x, y, time, conc, time_dbl)
    
    cca_prior_long <- bind_rows(cca_prior_long, pred_df_cca)
  }
  
  rccr_pred_long <- cca_prior_long |> 
    mutate(pred = "rCCR") |> 
    select(-time_dbl)
  
  rccr_pred_df <- rccr_pred_long |> 
    select(-pred) |> 
    pivot_wider(names_from = time, values_from = conc)
  
  list("df" = rccr_pred_df, "long" = rccr_pred_long)
}

# simple bayes
sb_pred_function <- function(well_xy, well_long){
  
  # well_xy = xy
  # well_long = long_well
  
  sb_df <- plume_xy
  
  for (tm in paste0("ts", 1:5)){ # initialize with the first five time steps
    well_time <- filter(well_long, time==tm)
    mod_tps <- Tps(x = select(well_time, x, y), 
                   Y = well_time$conc)
    pred <- predict(mod_tps, x = plume_xy)[,1]
    pred[pred < 0] <- 0
    sb_df <- sb_df |> 
      mutate(!!sym(tm) := pred)
  }
  
  sb_long_df <- sb_df |> 
    pivot_longer(starts_with("ts"), names_to = "time", values_to = "conc") |> 
    mutate(time_dbl = str_extract(time, "\\d+") |> as.numeric()) |> 
    filter(time_dbl <= 5) |> 
    left_join(mutate(well_xy, target = TRUE), by = c("x", "y")) |> 
    mutate(target = ifelse(is.na(target), FALSE, TRUE)) |> # locations with NA are not observation locations
    arrange(-target, x, y) # arranges in descending order, for boolean starts with TRUE
  
  d_mat <- sb_long_df |> 
    filter(time_dbl == 1) |> 
    arrange(-target, x, y) |> 
    select(x, y) |> 
    dist(upper = TRUE, diag = TRUE) |> 
    as.matrix()
  
  # get covariance matrix using Wendlend Covariance function
  phi <- Wendland(d_mat, k = 1, dimension = 2) |> as.matrix() 
  
  #for loop starts here
  for (i in 6:40){

    # get v and lambda at each grid location based on previous 5 time steps
    prev_stats <- sb_long_df |> 
      filter(time_dbl < i) |> 
      group_by(x, y, target) |> 
      summarize(v = mean(conc), lambda = sd(conc)) |> 
      ungroup() |> 
      arrange(-target, x, y)
    
    lambda <- diag(prev_stats$lambda) # matrix of standard deviations
    psi <- lambda %*% phi %*% lambda # scale matrix
    
    v1 <- prev_stats |> filter(target) |> pull(v)
    v2 <- prev_stats |> filter(!target) |> pull(v)
    psi2_1 <- psi[(length(v1) + 1):nrow(psi), 1:length(v1)]
    psi1_1 <- psi[1:length(v1), 1:length(v1)] +  matrix(
      rnorm(length(v1)^2, 0, sd = 0.0000001), 
      nrow = length(v1))
    
    x1 <- well_long |> 
      filter(time_dbl == i) |> 
      arrange(x, y) |> 
      pull(conc)
    
    v2_1 <- v2 + psi2_1 %*% solve(psi1_1, tol = 1e-17) %*% (x1 - v1)
    
    pred_surf <- tibble(
      x = prev_stats$x,
      y = prev_stats$y,
      conc = c(x1, v2_1),
      time = paste0("ts", i),
      time_dbl = i,
      target = prev_stats$target
    ) |> 
      mutate(conc = ifelse(conc < 0, 0, conc))
    
    ggplot(pred_surf, aes(x, y, fill = conc)) + geom_raster()
    
    sb_long_df <- bind_rows(sb_long_df, pred_surf)
  }
  
  sb_pred_long <- sb_long_df |> 
    mutate(pred = "SB") |> 
    select(x, y, time, conc, pred)
  
  sb_pred_df <- sb_pred_long |> 
    select(-pred) |> 
    pivot_wider(names_from = time, values_from = conc)
  
  list("df" = sb_pred_df, "long" = sb_pred_long)
}

# random forest -- one step at a time
rfsi_pred_function <- function(well_xy, well_long){
  
  rfsi_pred_df <- plume_xy
  tm_interest <- paste0("ts", 1:40)
  well_sub <- well_long |> dplyr::select(x, y, conc, time_dbl, time)
  
  # set station id for all coordinates
  plume_xy$id <- row_number(plume_xy)
  # join with plume_xy to get correct well id for model fitting
  well_sub <- well_sub |>
    right_join(plume_xy, by = c("x", "y"))

  for (tm in tm_interest){
    
    well_time <- filter(well_sub, time==tm)

    data = well_time
    data <- data |> dplyr::select(id, x, y, conc)
    data.staid.x.y.z <- c(1,2,3,NA) 
    class(data) <- "data.frame"
    
    fm.RFSI <- as.formula("conc ~ 1") # conc ~ 1 allows only nearest observations and distances to them
    
    
    mod_rfsi <- meteo::rfsi(formula = fm.RFSI,
                            data = data,
                            data.staid.x.y.z = data.staid.x.y.z,
                            n.obs = 5, # number of nearest observations 
                            cpus = detectCores()-1,
                            progress = TRUE,
                            importance = "impurity",
                            seed = 42,
                            num.trees = 250,
                            mtry = 5,
                            splitrule = "variance",
                            min.node.size = 5,
                            sample.fraction = 0.95,
                            quantreg = FALSE
                            )
    newdata <- data.frame(plume_xy)
    newdata <- newdata |> select(id, x, y) # rearrange to have id as first column
    newdata.staid.x.y.z <- c("id","x","y",NA)
    
    pred <- meteo::pred.rfsi(model = mod_rfsi,
                    data = data, 
                    data.staid.x.y.z = data.staid.x.y.z,
                    obs.col = "conc",
                    output.format = "data.frame",
                    newdata = newdata,
                    newdata.staid.x.y.z = newdata.staid.x.y.z,
                    zero.tol = 0,
                    cpus = 1, # detectCores()-1,
                    progress = TRUE) |>
      dplyr::pull(pred)
    
      # pred |> dplyr::filter(id == 316)
      
    rfsi_pred_df <- rfsi_pred_df |>
      dplyr::mutate(!!sym(tm) := pred)
    
  }
  
  list(
    "df" = rfsi_pred_df,
    "long" = rfsi_pred_df |> 
      pivot_longer(starts_with("ts"), names_to = "time", values_to = "conc") |> 
      mutate(pred = "rsfi")
  )
  
}



# INLA
inla_pred_function <- function(well_xy, well_long){
  
  # create sf object from prediction (regular) grid
  sf_pred_grid <- st_as_sf(plume_xy, coords = c("x", "y"))
  
  # To get the corner points as a new sf object:
  bbox_points <- st_as_sfc(st_bbox(sf_pred_grid)) %>%
    st_cast("POLYGON") %>%
    st_boundary() %>%
    st_cast("POINT")
  
  # get bounding box of prediction grid
  bbox_pred_grid <- st_bbox(sf_pred_grid) # function will apply to a bbox
  # print(bbox_pred_grid)
  mesh_coords <- generate_well_spaced_points(bbox_pred_grid, 200, 5)
  
  # get interpolated boundary points
  bbox_poly <- st_as_sfc(bbox_pred_grid, crs = "ESPG:3421")
  num_points_per_edge <- 10 # Adjust as needed
  
  interpolated_points <- interpolate_bbox(bbox_poly, num_points_per_edge)
  boundary_points <- rbind(interpolated_points, st_coordinates(bbox_points))
  
  # create mesh using mesh coordinates
  mesh = fmesher::fm_mesh_2d_inla(
    loc = mesh_coords,
    loc.domain = boundary_points, # points along the border of the domain
    max.edge = diff(range(plume_xy$x))*(1/3) # make smaller to 1/3
  )
  # plot(mesh)
  # points(well_xy)
  
  spde = inla.spde2.matern(mesh=mesh)
  
  n_time <- length(unique(well_long$time_dbl)) # total number of time points
  time_index <- as.numeric(factor(well_long$time_dbl)) # years as factors
  # need to create a combined index for the space-time field
  well_long$st_index <- 1:nrow(well_long)
  
  A.est = inla.spde.make.A(
    mesh = mesh,
    loc = as.matrix(well_long[,c("x","y")]), # observation coordinates
    n.group = n_time, # number of year options
    group = time_index # maps each row to associated year
    
  )
  
  field.indices = inla.spde.make.index(
    "field",
    n.spde = mesh$n,
    n.group = n_time
  )
  
  stack.est = inla.stack(
    data = list(conc = well_long$conc),
    A = list(A.est, 1),
    effects = list(spatial.field = field.indices, 
                   Intercept=rep(1, nrow(well_long))),
    tag = "est")
  
  formula <- (conc ~ -1 + Intercept + 
                f(field, 
                  model=spde,
                  group=field.group, 
                  control.group=list(model="ar1")
                )
  )
  
  # initial model optimization over parameters
  tic("INLA fit")
  mod.mode = inla(formula,
                  data=inla.stack.data(stack.est, spde=spde),
                  family="gaussian",
                  control.predictor=list(A=inla.stack.A(stack.est), compute=FALSE),
                  # verbose = FALSE
  )
  toc()
  # INLA fit: 34.959 sec elapsed
  
  conc_reg_grid_results <- matrix(NA, nrow = nrow(plume_xy), ncol = n_time)
  conc_obs_grid_results <- matrix(NA, nrow = nrow(well_long), ncol = n_time)
  
  # prediction for each year
  for(i_time in 1:n_time){ # i_time = 1
    # for(i_day in 1:45){ # for testing - use previous line after testing
    
    cat("INLA Prediction Time:", i_time, fill = TRUE)
    
    A.pred = inla.spde.make.A(
      mesh = mesh, 
      loc = as.matrix(plume_xy[,c("x","y")]),
      group = rep(i_time, # i_time is the prediction time step
                  nrow(plume_xy)), 
      n.group = n_time
    )
    
    stack.pred = inla.stack(
      data = list(conc=NA),
      A=list(A.pred,1),
      effects = list(spatial.field = field.indices, 
                     Intercept=rep(1, nrow(plume_xy))),
      tag="pred")
    
    stack = inla.stack(stack.est, stack.pred)
    
    
    # prediction at prediction points
    tic("INLA prediction")
    mod = inla(formula,
               data=inla.stack.data(stack, spde=spde),
               family="gaussian",
               control.predictor=list(A=inla.stack.A(stack), compute=TRUE),
               control.mode=list(theta=mod.mode$mode$theta, restart=FALSE)
    )
    toc()
    # INLA prediction: 4.395 sec elapsed (one time step)
    # summary(mod)
    
    # store results for prediction grid and estimation grid
    index.pred = inla.stack.index(stack,"pred")$data
    index.est = inla.stack.index(stack, "est")$data
    conc_reg_grid_results[,i_time] <- mod$summary.linear.predictor[index.pred, "mean"]
    conc_obs_grid_results[,i_time] <- mod$summary.linear.predictor[index.est, "mean"]
    
  }
  
  # set anything less than zero to zero
  conc_reg_grid_results[conc_reg_grid_results < 0] <- 0
  
  inla_pred_df <- cbind(plume_xy, conc_reg_grid_results)
  colnames(inla_pred_df) <- c("x", "y", paste0("ts", 1:n_time))
  
  # ggplot() +
  #   geom_tile(data = inla_pred_df, aes(x = x, y = y, fill = ts1)) +
  #   viridis::scale_fill_viridis(discrete = FALSE) +
  #   geom_point(data = well_xy, aes(x = x, y = y), col = "magenta")
  # 
  # output two matrices
  list(
    "df" = inla_pred_df,
    "long" = inla_pred_df |> 
      pivot_longer(starts_with("ts"), names_to = "time", values_to = "conc") |> 
      mutate(pred = "inla")
  )
  
  
}

# simulation functions ---------------------------------------------------------

# calculate error for each method -- predicts using methods above, returns error
model_mspe <- function(xy, long_well,
                       n_wells, plume_type,
                       network_id
                       ){
  
  # method predictions
  tic()
  tpm <- tpm_pred_function(well_xy = xy, well_long = long_well)
  tps <- tps_pred_function(xy, long_well)
  tps_time <- tps_time_pred_function(xy, long_well)
  rccr <- rccr_pred_function(xy, long_well)
  sb <- sb_pred_function(xy, long_well)
  rfsi <- rfsi_pred_function(xy, long_well)
  inla <- inla_pred_function(well_xy = xy, well_long = long_well)
  toc()
  
  tpm_na <- tpm$long |> 
    filter(time == "ts10") |>
    mutate(tpm_na = is.na(conc)) |> 
    select(x, y, tpm_na)
  
  pred_df_long <- bind_rows(tps$long, tps_time$long, 
                            tpm$long, rccr$long,
                            sb$long, rfsi$long,
                            # uk$long, ukt$long,
                            inla$long,
                            # ebk$long # excluding this for now because it's taking forever and need to test
                            ) |> 
    #filter(time %in% tm_interest) |> 
    left_join(
      select(plume_data_long, x, y, time, obs_conc = conc),
      by = c("x", "y", "time")
    ) |> 
    mutate(
      conc = ifelse(conc < 0, 0, conc), #natural boundary of 0
      spe  = (conc - obs_conc)^2
    )
  
  full_grid <- pred_df_long |> 
    filter(pred != "tpm") |> 
    group_by(time, pred) |> 
    summarize(mspe = mean(spe, na.rm = TRUE)) |> 
    ungroup() |> 
    mutate(full_grid = TRUE)
  
  tpm_grid <- pred_df_long |> 
    left_join(tpm_na, by = c("x", "y")) |> 
    filter(!tpm_na) |> 
    group_by(time, pred) |> 
    summarize(mspe = mean(spe, na.rm=TRUE)) |> 
    ungroup() |> 
    mutate(full_grid = FALSE)
  

  all_mspe <- bind_rows(full_grid, tpm_grid)
  
 
  # prediction values for TPM without NAs
  tpm_pred_no_na <- tpm$long |>
    drop_na() |>
    left_join(
      select(plume_data_long, x, y, time, obs_conc = conc),
      by = c("x", "y", "time")
    )
  # combine with all other methods
 pred_df_out <-  pred_df_long |> 
    select(x, y, time, conc, pred, obs_conc) |>
    bind_rows(tpm_pred_no_na)

  
  # save prediction results -- all with full grid and tpm within convex hull
  save(pred_df_out,
       file = file.path("simulation","data", "results",
                        paste0(plume_type, "_pred_plume_network_", network_id, "_", n_wells, "_nwells.RDS"))
  )
  
  return(all_mspe)
  
}


# function to get well network correlations to zero
check_cond <- function(well_xy, well_long){
  
  cca_prior <- plume_xy
  
  for (i in 1:5){
    mod_tps <- Tps(x = well_xy,
                   Y = pull(filter(well_long, time_dbl == i), conc),
                   verbose = FALSE)
    pred <- predict(mod_tps, x = plume_xy)[,1]
    pred[pred < 0] <- 0
    cca_prior <- cca_prior |>
      mutate(!!sym(paste0("ts",i)) := pred)
  }
  
  cca_prior_long <- cca_prior |>
    pivot_longer(starts_with("ts"), names_to = "time", values_to = "conc") |>
    mutate(time_dbl = str_extract(time, "\\d+") |> as.numeric())
  
  x_train <- well_long |>
    select(x, y, conc, time_dbl) |>
    filter(time_dbl < 6) |>
    pivot_wider(names_from = c("x", "y"), values_from = "conc") |>
    select(-time_dbl) |>
    as.matrix()
  
  
  all_train <- cca_prior_long |>
    select(-time) |>
    filter(time_dbl < 6) |>
    pivot_wider(names_from = c("x", "y"), values_from = "conc")  |>
    select(-time_dbl) |>
    as.matrix()
  
  all_train <- all_train + matrix(
    rnorm(length(all_train), 0, sd = 0.000000000001), 
    nrow = nrow(all_train))
  
  cca_mod_all <- scca(x_train, all_train, scale = TRUE, verbose = FALSE)
}

# function to simulate network of observation wells
simulate_well_network <- function(i, well_n_int, plume_data_wide, plume_type, plume_data_long){
  
  well_v <- sort(well_n_int, decreasing = TRUE)
  
  well_n <- well_v[1] # start with first number of wells in network

  all_good <- FALSE
  
  # simulate well network
  while (!all_good) {

      # get starting position of plume
      plume_xy_start <- plume_data_wide |>
        filter(ts1 > 1) |>
        select(x, y)
  
      # get ending position of plume and simulate network so that
      #   65% of wells are in the plume area
      plume_xy_end <- plume_data_wide |>
        filter(ts40 > 1) |>
        select(x, y)
      
      close_min_y <- min(plume_xy_end$y)
      close_max_y <- max(plume_xy_end$y)
      
    n_close <- round(0.65 * well_n) # select number of close wells
    n_far <- well_n - n_close # remaining number of far wells
    
    well_obs_close <- plume_data_wide |>
      filter(y >= close_min_y & y <= close_max_y) |> # set y-range for close wells
      slice_sample(n = n_close) # randomly select rows

    well_obs_far <- plume_data_wide |>
      filter(y > close_max_y | y < close_min_y) |>
      slice_sample(n = n_far) # randomly select far wells
    
    well_obs <- bind_rows(well_obs_close, well_obs_far)
    
    well_xy <- select(well_obs, x, y)
    
    # keep simulated data only at chosen well locations:
    well_long <- plume_data_long |> 
      inner_join(well_xy)
    
    all_good <- !berryFunctions::is.error(check_cond(well_xy, well_long))
  }
  
  # get results for first set of well networks
  result_df <- model_mspe(
    well_xy, 
    well_long,
    plume_type = plume_type,
    network_id = i, 
    n_wells = well_n 
    ) |> 
    mutate(n_wells = well_n)
  

  # for the other number of wells, remove number of wells from max total
  for (wn in well_v[-1]){
    print(wn)
    all_good <- FALSE
    while (!all_good) {
      n_close <- round(0.65 * wn)
      n_far <- wn - n_close
      
      well_obs_close <- well_obs |>
        filter(y >= close_min_y & y <= close_max_y) |> # set y-range for close wells
        slice_sample(n = n_close) # randomly select rows

      well_obs_far <- well_obs |>
        filter(y > close_max_y | y < close_min_y) |>
        slice_sample(n = n_far) # randomly select far wells

      well_n <- bind_rows(well_obs_close, well_obs_far)

      well_xy_n <- select(well_n, x, y)

      well_long_n <- plume_data_long |>
        inner_join(well_xy_n, by = c("x", "y"))

      all_good <- !berryFunctions::is.error(check_cond(well_xy_n, well_long_n))
    }
    well_obs <- well_n
    well_xy <- well_xy_n
    well_long <- well_long_n

    result_df <- model_mspe( well_xy, 
                             well_long,
                             plume_type = plume_type,
                             network_id = i, 
                             n_wells = wn 
                             ) |>
      mutate(n_wells = wn) |>
      bind_rows(result_df)
  }

  # ID this well network with i
  result_df |> mutate(network_id = i)
  
}
