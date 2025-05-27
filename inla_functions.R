library(INLA)
library(sf)

# INLA functions

generate_adaptive_mesh <- function(data_points, high_density_threshold,
                                   min_edge_high, max_edge_high,
                                   max_edge_low, cutoff_low,
                                   boundary = NULL) {
  if (st_geometry_type(data_points)[1] != "POINT") {
    stop("data_points must be an sf object of POINT geometries.")
  }
  
  # Get the CRS of the data_points object
  data_crs <- st_crs(data_points)
  
  # identify high-density areas
  high_density_coords <- st_coordinates(data_points)
  is_high_density <- apply(high_density_coords, 1, function(pt) {
    current_point_sfc <- st_sfc(st_point(pt), crs = data_crs)
    distances <- sf::st_distance(current_point_sfc, data_points)
    distance_values <- as.numeric(distances)
    threshold <- as.numeric(mean(st_bbox(data_points)[c("xmax", "ymax")] -
                                   st_bbox(data_points)[c("xmin", "ymin")])) / 20
    return(sum(distance_values <= threshold) > high_density_threshold)
  })
  
  print("Head of is_high_density:")
  print(head(is_high_density))
  print(paste("Length of is_high_density:", length(is_high_density)))
  print(paste("Number of rows in data_points:", nrow(data_points)))
  
  high_density_subset <- data_points[is_high_density, ]
  low_density_subset <- data_points[!is_high_density, ]
  
  print(paste("Number of high-density points:", nrow(high_density_subset)))
  print(paste("Number of low-density points:", nrow(low_density_subset)))
  
  # generate interior knots with varying density
  n_high <- nrow(high_density_subset)
  n_low <- nrow(low_density_subset)
  
  interior_knots_high <- NULL
  if (n_high > 0) {
    bbox_high <- st_bbox(high_density_subset)
    ik_high_df <- generate_well_spaced_points(
      bbox_high[c("xmin", "ymin", "xmax", "ymax")],
      n_high * 2,
      min_edge_high / 2
    )
    interior_knots_high <- st_as_sf(ik_high_df, coords = c("x", "y"), crs = data_crs)
  }
  
  interior_knots_low <- NULL
  if (n_low > 0) {
    bbox_low <- st_bbox(low_density_subset)
    ik_low_df <- generate_well_spaced_points(
      bbox_low[c("xmin", "ymin", "xmax", "ymax")],
      n_low * 1.5,
      max_edge_low / 2
    )
    interior_knots_low <- st_as_sf(ik_low_df, coords = c("x", "y"), crs = data_crs)
  }
  
  all_interior_knots_sf <- rbind(interior_knots_high, interior_knots_low)
  
  if (is.null(all_interior_knots_sf)) {
    bbox_overall <- st_bbox(data_points)
    ik_overall_df <- generate_well_spaced_points(
      bbox_overall[c("xmin", "ymin", "xmax", "ymax")],
      nrow(data_points) * 0.5,
      max_edge_low / 2
    )
    all_interior_knots_sf <- st_as_sf(ik_overall_df, coords = c("x", "y"), crs = data_crs)
  }
  
  # 3. Create the adaptive mesh with varying max.edge
  adaptive_mesh <- inla.mesh.2d(
    loc = st_coordinates(all_interior_knots_sf),
    boundary = boundary,
    max.edge = c(max_edge_high, max_edge_low),
    cutoff = cutoff_low
  )
  
  return(adaptive_mesh)
}

# helper function for well-spaced points
generate_well_spaced_points <- function(bbox, n_points, min_dist) {
  points <- data.frame(x = runif(1, bbox[1], bbox[3]),
                       y = runif(1, bbox[2], bbox[4]))
  
  while (nrow(points) < n_points * 2) {
    new_point <- data.frame(x = runif(1, bbox[1], bbox[3]),
                            y = runif(1, bbox[2], bbox[4]))
    
    dists <- spDistsN1(as.matrix(points[, c("x", "y")]), as.numeric(new_point))
    if (all(dists > min_dist)) {
      points <- rbind(points, new_point)
    }
    if (nrow(points) >= n_points) break
  }
  
  if (nrow(points) > n_points) {
    points <- points[sample(nrow(points), n_points), ]
  }
  
  return(points)
}



interpolate_bbox <- function(poly, n) {
  boundary_lines <- st_cast(poly, "LINESTRING")
  coords <- st_coordinates(boundary_lines)
  num_segments <- nrow(coords) - 1
  interpolated_coords_list <- list()
  
  for (i in 1:num_segments) {
    start <- coords[i, 1:2]
    end <- coords[i + 1, 1:2]
    segment <- st_linestring(matrix(c(start, end), nrow = 2, byrow = TRUE))
    interpolated <- st_line_sample(segment, n = n)
    interpolated_coords <- st_coordinates(interpolated)[, 1:2]
    interpolated_coords_list <- append(interpolated_coords_list, list(interpolated_coords))
  }
  return(do.call(rbind, interpolated_coords_list))
}
