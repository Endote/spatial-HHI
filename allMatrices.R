# --- Libraries ---
libs <- c("sf", "dplyr")
invisible(lapply(libs, function(pkg) if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)))
lapply(libs, library, character.only = TRUE)

# --- Manual grid locations in [0,1] space (your original recipe) ---
manual_locations <- list(
  n2 = data.frame(
    x = c(0.25, 0.75),
    y = c(0.25, 0.75)
  ),
  n4 = data.frame(
    x = c(0.25, 0.25, 0.75, 0.75),
    y = c(0.25, 0.75, 0.25, 0.75)
  ),
  n8 = expand.grid(
    x = seq(0.125, 0.875, length.out = 4),
    y = c(0.25, 0.75)
  ),
  n16 = expand.grid(
    x = seq(0.125, 0.875, length.out = 4),
    y = seq(0.125, 0.875, length.out = 4)
  )
)
manual_locations$n16 <- manual_locations$n16[order(manual_locations$n16$y, manual_locations$n16$x), ]

# --- Precompute all random locations (so cols 3/4 are identical per n) ---
random_school_locations <- list()
for (n in c(2,4,8,16)) {
  set.seed(222)
  random_school_locations[[paste0("n", n)]] <- data.frame(
    x = runif(n, min=0, max=1),
    y = runif(n, min=0, max=1)
  )
}

# --- Params for all plots ---
n_pop_side <- 32
total_pop <- 5000
boundary <- st_polygon(list(rbind(c(0,0), c(1,0), c(1,1), c(0,1), c(0,0)))) %>% st_sfc(crs = 4326)
crs_proj <- 3857
boundary_proj <- st_transform(boundary, crs_proj)

# Only show n=2, 4, 8
row_n <- c(2, 4, 8)
col_scenarios <- c("Uniform Loc + Uniform Pop", "Uniform Loc + Random Pop", "Random Loc + Uniform Pop", "Random Loc + Random Pop")

# --- Main plot function ---
plot_voronoi_scenario <- function(n, loc_type, pop_type, manual_locations, random_school_locations, irow, icol) {
  # -- School points --
  if (loc_type == "uniform") {
    coords <- manual_locations[[paste0("n", n)]]
    school_pts <- st_as_sf(coords, coords = c("x", "y"), crs = 4326)
    school_pts_proj <- st_transform(school_pts, crs_proj)
  } else {
    coords <- random_school_locations[[paste0("n", n)]]
    school_pts <- st_as_sf(coords, coords = c("x", "y"), crs = 4326)
    school_pts_proj <- st_transform(school_pts, crs_proj)
  }
  
  # -- Population tiles (always a grid) --
  pop_tiles <- st_make_grid(boundary, n = c(n_pop_side, n_pop_side)) %>% st_sf()
  pop_tiles_proj <- st_transform(pop_tiles, crs_proj)
  n_tiles <- nrow(pop_tiles_proj)
  
  # -- Population assignment --
  if (pop_type == "uniform") {
    pop_per_tile <- total_pop / n_tiles
    pop_tiles_proj$pop <- rep(pop_per_tile, n_tiles)
  } else {
    set.seed(200 + irow*10 + icol) # Reproducibility
    base_pop <- 1
    lambda <- 10
    raw_pop <- base_pop + rpois(n_tiles, lambda = lambda)
    n_clusters <- min(4, n)
    cluster_centers <- sample(1:n_tiles, n_clusters)
    for(center in cluster_centers) {
      dists <- st_distance(pop_tiles_proj[center, ], pop_tiles_proj)
      dists_num <- as.numeric(dists)
      cluster_boost <- as.numeric(80 * exp(-((dists_num/150)^2)))
      raw_pop <- raw_pop + round(cluster_boost)
    }
    pop_tiles_proj$pop <- round(raw_pop / sum(raw_pop) * total_pop)
    pop_tiles_proj$pop[pop_tiles_proj$pop == 0] <- 1
    pop_tiles_proj$pop <- round(pop_tiles_proj$pop / sum(pop_tiles_proj$pop) * total_pop)
    diff <- total_pop - sum(pop_tiles_proj$pop)
    if (diff != 0) pop_tiles_proj$pop[1] <- pop_tiles_proj$pop[1] + diff
  }
  pop_tiles_proj$tile_id <- 1:nrow(pop_tiles_proj)
  
  # -- Voronoi tessellation --
  vor <- st_voronoi(st_union(school_pts_proj), envelope = boundary_proj)
  vor_polys <- st_collection_extract(vor, "POLYGON")
  vor_sf <- st_sf(geometry = vor_polys)
  vor_sf <- st_intersection(vor_sf, boundary_proj)
  vor_sf$school_id <- 1:nrow(vor_sf)
  
  # -- Fractional assignment (best for accuracy) --
  pop_vor <- st_intersection(pop_tiles_proj, vor_sf)
  pop_vor$sub_area <- as.numeric(st_area(pop_vor))
  tile_areas <- as.numeric(st_area(pop_tiles_proj))
  names(tile_areas) <- pop_tiles_proj$tile_id
  pop_vor$tile_area <- tile_areas[as.character(pop_vor$tile_id)]
  pop_vor$pop_share <- pop_vor$pop * (pop_vor$sub_area / pop_vor$tile_area)
  
  pop_sum <- pop_vor %>%
    st_drop_geometry() %>%
    group_by(school_id) %>%
    summarize(pop_tile = sum(pop_share, na.rm=TRUE), .groups = 'drop')
  vor_sf <- left_join(vor_sf, pop_sum, by = "school_id")
  vor_sf$pop_tile[is.na(vor_sf$pop_tile)] <- 0
  
  # -- Spatial HHI calculation --
  Pi <- vor_sf$pop_tile
  Ai <- as.numeric(st_area(vor_sf))
  Pi[is.na(Pi)] <- 0
  Ai[is.na(Ai)] <- 1e-4
  pi <- Pi / Ai
  wi <- pi / sum(pi, na.rm=TRUE)
  Ptot <- sum(Pi, na.rm=TRUE)
  ms <- Pi / Ptot
  S_i <- ms * wi
  S_i <- S_i / sum(S_i, na.rm=TRUE)
  vor_sf$spatially_weighted_share <- S_i
  hhi_val <- sum(S_i^2, na.rm=TRUE)
  
  # -- ROUNDED COLOR BINS LOGIC --
  share_vals <- round(vor_sf$spatially_weighted_share, 4)
  unique_shares <- sort(unique(share_vals))
  pal <- colorRampPalette(c("yellow", "purple"))(length(unique_shares))
  vor_sf$color <- pal[match(share_vals, unique_shares)]
  
  # -- PLOT: Each subplot
  plot(
    st_geometry(vor_sf), col = as.character(vor_sf$color), border="gray30",
    main = NA, axes = FALSE
  )
  plot(st_geometry(school_pts_proj), col="black", bg="green", pch=21, cex=1.2, add=TRUE)
  
  # --- LEGEND for every plot, bigger and to the right ---
  legend_labels <- sprintf("%.3f", unique_shares)
  legend(
    x = "topright", 
    inset = c(0.8, 0),  # Further right, tweak as needed
    legend = legend_labels,
    fill = pal,
    border = "gray30",
    title = "Pop Share",
    bty = "n",
    cex = 1.15,
    xpd = TRUE,
    horiz = FALSE,
  )
  
  # --- Annotate Spatial HHI and title ---
  title(
    main = sprintf("%s\nn = %d\nSpatial HHI: %.3f", 
                   col_scenarios[icol], n, hhi_val),
    line = -0.5, cex.main = 0.9
  )
}

# ---- Plotting 3x4 grid (no 4th row) ----
par(mfrow = c(3,4), mar = c(1.5, 2.5, 3.8, 1.5), oma = c(3, 3, 3, 1))
for (i in 1:3) {
  n <- row_n[i]
  for (j in 1:4) {
    if (j == 1) plot_voronoi_scenario(n, "uniform", "uniform", manual_locations, random_school_locations, i, j)
    if (j == 2) plot_voronoi_scenario(n, "uniform", "random", manual_locations, random_school_locations, i, j)
    if (j == 3) plot_voronoi_scenario(n, "random",  "uniform", manual_locations, random_school_locations, i, j)
    if (j == 4) plot_voronoi_scenario(n, "random",  "random",  manual_locations, random_school_locations, i, j)
  }
}

mtext("Behaviour of Spatially Weighted HHI, Number of Points (2,4,8)", side = 3, line = 0.5, outer = TRUE, cex = 1.6, font = 2)
