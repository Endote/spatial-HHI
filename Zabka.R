# ---- SETUP ----

# Install if missing
pkgs <- c("osmdata", "dplyr", "sf", "RColorBrewer", "spatstat", "spatstat.geom", "tidyr", "spdep", "igraph", "GenSA", "future.apply")
invisible(lapply(pkgs, function(pkg) if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)))
lapply(pkgs, library, character.only = TRUE)

options(scipen=999)

plan(multisession)

compute_spatially_weighted_hhi <- function(vor_sf) {
  # Pi = population in tile (RES)
  Pi <- vor_sf$pop_tile
  # Ai = area of tile, in km²
  Ai <- as.numeric(st_area(vor_sf)) / 1e6
  Pi[is.na(Pi)] <- 0
  Ai[is.na(Ai)] <- 1e-4 # avoid div by 0
  
  # pi = population density (persons/km²)
  pi <- Pi / Ai
  
  # wi = normalized population density
  wi <- pi / sum(pi, na.rm=TRUE)
  
  # Market share by population
  Ptot <- sum(Pi, na.rm=TRUE)
  ms <- Pi / Ptot
  
  # Spatially weighted share (product, then normalize)
  S_i <- ms * wi
  S_i <- S_i / sum(S_i, na.rm=TRUE)
  
  vor_sf$spatially_weighted_share <- S_i
  
  # HHI = sum of squares of spatially weighted shares
  hhi_sw <- sum(S_i^2, na.rm=TRUE)
  
  list(
    vor_sf = vor_sf,
    hhi_spatial_weighted = hhi_sw
  )
}

# ---- DATA LOAD ----

setwd('/Users/norbert.jaworski/Documents/uni/SpatialR') # Change as needed

POW <- st_read("data/powiaty/powiaty.shp", stringsAsFactors = FALSE)
POW$jpt_nazwa_ <- iconv(POW$jpt_nazwa_, "latin1", "UTF-8")
POW <- st_transform(POW, 4326)
POW.waw <- POW[POW$jpt_nazwa_ == "powiat Warszawa",]

WAW <- st_read("data/dzielnice_Warszawy/dzielnice_Warszawy.shp") %>% st_transform(4326)
POP.WAW <- st_read("data/POP.WAW/POP.WAW.shp")

warsaw_bbox <- getbb("Warsaw, Poland")

river <- warsaw_bbox %>%       
  opq() %>%
  add_osm_feature(key="waterway", value="river", bbox = warsaw_bbox) %>%
  osmdata_sf()

railway <- warsaw_bbox %>%       
  opq() %>%
  add_osm_feature(key="railway", value='rail', bbox = warsaw_bbox) %>%
  osmdata_sf()

# ---- 1. SCHOOL POINTS: GET & CLEAN ----

get_zabka_points <- function(bbox, waw_sf, threshold_m = 200) {
  message("Fetching OSM school points...")
  
  zabkas <- osmdata_sf(
    opq(bbox = bbox) %>%
      add_osm_feature(key = "shop") %>% 
    add_osm_feature(key = "name", value = "Żabka")
  )$osm_points
  
  
  # Defensive: Check if any points returned
  if (is.null(zabkas) || nrow(zabkas) == 0) {
    warning("No Żabka points returned from OSM in this bounding box.")
    return(NULL)
  }
  
  zabka_points <- st_intersection(zabkas, waw_sf)
  
  if (nrow(zabka_points) == 0) {
    warning("No Żabka shops within the Warsaw boundary.")
    return(NULL)
  }
  
  # Medoid-based deduplication
  merge_points_medoid <- function(points_sf, threshold_m = 100, metric_crs = 2180, output_crs = 4326) {
    if (nrow(points_sf) == 0) return(points_sf)
    points_metric <- st_transform(points_sf, metric_crs)
    dist_matrix <- st_distance(points_metric)
    threshold <- units::set_units(threshold_m, "m")
    adj_matrix <- (dist_matrix < threshold) & (dist_matrix > units::set_units(0, "m"))
    g <- igraph::graph_from_adjacency_matrix(as.matrix(adj_matrix), mode = "undirected", diag = FALSE)
    clusters <- igraph::components(g)$membership
    points_metric$cluster <- clusters
    medoid_indices <- c()
    for (cl in sort(unique(clusters))) {
      idx <- which(points_metric$cluster == cl)
      if (length(idx) == 1) {
        medoid_indices <- c(medoid_indices, idx)
      } else {
        sub_dm <- as.matrix(dist_matrix[idx, idx])
        sum_dists <- rowSums(sub_dm)
        medoid_local_index <- which.min(sum_dists)
        medoid_indices <- c(medoid_indices, idx[medoid_local_index])
      }
    }
    medoid_points <- points_metric[unique(medoid_indices), ]
    st_transform(medoid_points, output_crs)
  }
  
  merge_points_medoid(zabka_points, threshold_m)
}

zabka_medoids <- get_zabka_points(warsaw_bbox, WAW, 200)


# ---- 2. VORONOI ----

get_voronoi <- function(points_sf, boundary_sf, crs_metric = 2180) {
  pts <- st_transform(points_sf, crs_metric)
  boundary <- st_transform(boundary_sf, crs_metric)
  envelope <- st_union(boundary)
  voronoi_geom <- pts %>% st_geometry() %>% st_union() %>% st_voronoi(envelope = st_geometry(envelope))
  voronoi_polygons <- st_collection_extract(voronoi_geom, "POLYGON") %>% st_sf(crs = st_crs(boundary))
  st_intersection(voronoi_polygons, boundary)
}

voronoi_clipped <- get_voronoi(zabka_medoids, POW.waw)

# ---- 3. POPULATION OVERLAY ----

get_tile_population <- function(voronoi_sf, pop_grid_sf, crs_metric = 2180) {
  POP_2180 <- st_transform(pop_grid_sf, crs_metric) %>%
    mutate(cell_id = row_number(), cell_area = st_area(geometry))
  VOR_2180 <- st_transform(voronoi_sf, crs_metric) %>%
    mutate(tile_id = row_number(), tile_area = st_area(geometry))
  
  pop_vor <- st_intersection(POP_2180, VOR_2180) %>%
    mutate(
      sub_area = st_area(geometry),
      pop_share = RES * (sub_area / cell_area)
    )
  pop_per_tile <- pop_vor %>%
    st_drop_geometry() %>%
    group_by(tile_id) %>%
    summarize(pop_tile = sum(pop_share, na.rm = TRUE))
  VOR_2180 <- VOR_2180 %>%
    left_join(pop_per_tile, by = "tile_id") %>%
    mutate(cell_area_km2 = as.numeric(tile_area) / 1e6,
           pop_density = as.numeric(pop_tile) / cell_area_km2)
  VOR_2180
}

VOR_2180 <- get_tile_population(voronoi_clipped, POP.WAW)

# ---- 4. SCHOOL ID ASSIGNMENT ----

assign_school_ids <- function(vor_sf, points_sf) {
  points_sf$school_id <- 1:nrow(points_sf)
  st_join(vor_sf, st_transform(points_sf %>% select(school_id), st_crs(vor_sf)), join = st_contains)
}
VOR_2180 <- assign_school_ids(VOR_2180, zabka_medoids)

# ---- 5. MARKET SHARE CALC ----

get_market_shares <- function(vor_sf) {
  total_pop <- sum(vor_sf$pop_tile, na.rm = TRUE)
  school_shares <- vor_sf %>%
    st_drop_geometry() %>%
    group_by(school_id) %>%
    summarize(pop_tile = sum(pop_tile, na.rm = TRUE)) %>%
    mutate(market_share = pop_tile / total_pop)
  vor_sf <- vor_sf %>% left_join(school_shares, by = "school_id")
  list(vor_sf = vor_sf, school_shares = school_shares)
}

market_out <- get_market_shares(VOR_2180)
VOR_2180 <- market_out$vor_sf
school_shares <- market_out$school_shares

# ---- 6. VISUALIZATION ----

plot_market_share_map <- function(vor_sf, school_pts, pow_sf) {
  VOR_4326 <- st_transform(vor_sf, 4326)
  school_pts_wgs <- st_transform(st_as_sf(school_pts), 4326)
  share_vals <- VOR_4326$spatially_weighted_share
  n_breaks <- 6
  share_breaks <- seq(min(share_vals, na.rm = TRUE), max(share_vals, na.rm = TRUE), length.out = n_breaks + 1)
  share_colors <- colorRampPalette(c("yellow", "purple"))(n_breaks)
  VOR_4326$color <- cut(VOR_4326$spatially_weighted_share, breaks = share_breaks, labels = share_colors, include.lowest = TRUE)
  plot(st_geometry(pow_sf), col = "white", border = "gray30", lwd = 1.5,
       main = "Zabka Spatially Weighted Market Share (Population Density)")
  plot(st_geometry(VOR_4326), col = as.character(VOR_4326$color), border = "black", lwd = 1.2, add = TRUE)
  plot(st_geometry(school_pts_wgs), col = "black", bg = "green", pch = 21, cex = 1.3, add = TRUE)
  plot(st_geometry(river$osm_lines), col = "steelblue1", lwd = 3, add = TRUE)
  plot(st_geometry(railway$osm_lines), col = "gray", lwd = 1, add = TRUE)
  legend("bottomleft",
         legend = paste0(round(share_breaks[-length(share_breaks)], 5), " to ", round(share_breaks[-1], 5)),
         fill = share_colors,
         title = "Spatially Weighted Share",
         bty = "n", cex = 0.8)
}

VOR_2180$pop_tile <- VOR_2180$pop_tile.x
sw_before <- compute_spatially_weighted_hhi(VOR_2180)
VOR_2180 <- sw_before$vor_sf

plot_market_share_map(VOR_2180, zabka_medoids, POW.waw)

# ---- 7. EXPANSION CANDIDATES ----

# --- PARAMETERS ---
n_top_tiles <- 10        # How many top market share tiles to scan
n_points_per_tile <- 6  # How many candidate points per tile

# --- 1. Get Top Market Share Tiles ---
top_tiles <- VOR_2180 %>%
  arrange(desc(market_share)) %>%
  slice_head(n = n_top_tiles)

# --- 2. Advanced Candidate Point Sampling ---
library(lhs)           # For Latin Hypercube

remove_near_border <- function(tile, pts, min_dist = 50) {
  dists <- st_distance(pts, st_boundary(tile))
  pts[dists > min_dist, , drop=FALSE]
}

# Latin Hypercube sampling within a polygon
sample_lhs_points <- function(poly, n) {
  bb <- st_bbox(poly)
  # Latin Hypercube in [0,1]^2 then mapped to bbox
  lhs_norm <- randomLHS(n, 2)
  xs <- bb["xmin"] + lhs_norm[,1] * (bb["xmax"] - bb["xmin"])
  ys <- bb["ymin"] + lhs_norm[,2] * (bb["ymax"] - bb["ymin"])
  pts_sf <- st_as_sf(data.frame(x=xs, y=ys), coords = c("x", "y"), crs = st_crs(poly))
  pts_sf[st_within(pts_sf, poly, sparse = FALSE), ]
}

sampling_method <- "lhs"  # options: "random", "lhs"
candidate_points <- do.call(rbind, lapply(1:nrow(top_tiles), function(i) {
  tile <- st_geometry(st_transform(top_tiles[i, ], 2180))[[1]]  # POLYGON only
  crs_tile <- st_crs(tile)
  if (sampling_method == "lhs") {
    candidates <- sample_lhs_points(tile, n_points_per_tile * 2)
  } else {
    candidates <- st_sample(tile, n_points_per_tile * 2, type = "random")
    if (is.null(st_crs(candidates))) st_crs(candidates) <- crs_tile
    candidates <- st_as_sf(st_sf(geometry = candidates)) # wrap as sf
  }
  candidates <- remove_near_border(tile, candidates, min_dist = 50)
  if (nrow(candidates) > n_points_per_tile) candidates <- candidates[1:n_points_per_tile,]
  if (nrow(candidates) > 0) {
    pts_sf <- st_sf(tile_id = top_tiles$tile_id[i], geometry = st_geometry(candidates))
    st_crs(pts_sf) <- crs_tile
    pts_sf$used_centroid <- FALSE
    pts_sf
  } else {
    centroid <- st_centroid(st_sfc(tile, crs = crs_tile))
    pts_sf <- st_sf(tile_id = top_tiles$tile_id[i], geometry = centroid)
    st_crs(pts_sf) <- crs_tile
    pts_sf$used_centroid <- TRUE
    pts_sf
  }
}))
candidate_points <- st_sf(candidate_points)
st_crs(candidate_points) <- 2180

# --- 3. Evaluate Each Candidate by Simulating School Addition ---

evaluate_candidate_impact <- function(candidate_point, zabka_pts, pop_grid, boundary_sf, crs_metric = 2180) {
  zabka_pts_metric <- st_transform(st_as_sf(zabka_pts), crs_metric)
  candidate_point_metric <- st_transform(candidate_point, crs_metric)
  zabka_coords <- st_coordinates(zabka_pts_metric)
  candidate_coords <- st_coordinates(candidate_point_metric)[1, 1:2]
  all_coords <- unique(rbind(zabka_coords, candidate_coords))
  multipoint <- st_multipoint(all_coords)
  multipoint_geom <- st_sfc(multipoint, crs = crs_metric)
  boundary <- st_union(st_transform(boundary_sf, crs_metric))
  vor_geom <- st_voronoi(multipoint_geom, envelope = st_geometry(boundary))
  vor_polys <- st_collection_extract(vor_geom, "POLYGON")
  if(length(vor_polys) != nrow(all_coords)){
    return(list(market_share = NA, voronoi = NULL))
  }
  vor_sf <- st_sf(geometry = vor_polys, crs = crs_metric)
  vor_clipped <- st_intersection(vor_sf, boundary)
  centroids <- st_centroid(vor_clipped)
  points_sf <- st_sf(zabka_id = 1:nrow(all_coords), geometry = st_sfc(lapply(1:nrow(all_coords), function(i) st_point(all_coords[i,])), crs = crs_metric))
  nearest_ids <- st_nearest_feature(centroids, points_sf)
  vor_clipped$zabka_id <- nearest_ids
  POP_2180 <- st_transform(pop_grid, crs_metric) %>%
    mutate(cell_id = row_number(), cell_area = st_area(geometry))
  vor_clipped <- vor_clipped %>% mutate(tile_id = row_number(), tile_area = st_area(geometry))
  pop_vor <- st_intersection(POP_2180, vor_clipped) %>%
    mutate(
      sub_area = st_area(geometry),
      pop_share = RES * (sub_area / cell_area) # (or RES, if full population)
    )
  pop_per_tile <- pop_vor %>%
    st_drop_geometry() %>%
    group_by(tile_id) %>%
    summarize(pop_tile = sum(pop_share, na.rm = TRUE))
  vor_clipped <- vor_clipped %>%
    left_join(pop_per_tile, by = "tile_id")
  
  # Compute the new candidate's total population tile
  # The candidate will be the *last* point (so max(zabka_id))
  candidate_id <- max(vor_clipped$zabka_id, na.rm=TRUE)
  total_pop <- sum(vor_clipped$pop_tile, na.rm = TRUE)
  candidate_pop <- vor_clipped$pop_tile[vor_clipped$zabka_id == candidate_id]
  market_share <- candidate_pop / total_pop
  
  return(list(
    market_share = market_share,
    voronoi = vor_clipped
  ))
}


# --- 4. Global HHI Optimization via Simulated Annealing ---

library(GenSA)

# Objective: negative of the new candidate's market share (since we minimize)
market_share_objective <- function(coord, zabka_pts, pop_grid, boundary_sf, city_poly) {
  pt <- st_sfc(st_point(coord), crs = 2180)
  poly_proj <- st_transform(city_poly, st_crs(pt))
  if (!as.logical(st_within(pt, poly_proj, sparse = FALSE)[1,1])) return(1e6)
  candidate_sf <- st_sf(geometry = pt)
  res <- evaluate_candidate_impact(candidate_sf, zabka_pts, pop_grid, boundary_sf)
  if (is.null(res$voronoi) || is.na(res$market_share)) return(1e6)
  # Return *negative* share to maximize share
  return(-as.numeric(res$market_share))
}

# Use all LHS candidates as seeds
city_poly <- st_union(POW.waw)
city_poly <- st_transform(city_poly, 2180)   # match everything

# Filter candidate_points to only those inside city
inside_city <- st_within(candidate_points, city_poly, sparse = FALSE)[,1]
cat("Candidate points inside city:", sum(inside_city), "/", nrow(candidate_points), "\n")
candidate_points <- candidate_points[inside_city,]
seed_coords <- st_coordinates(candidate_points)

best_hhi <- Inf
best_coord <- NULL
best_result <- NULL

bb <- st_bbox(city_poly)
lower <- c(bb["xmin"], bb["ymin"])
upper <- c(bb["xmax"], bb["ymax"])

# Prepare a function for one seed
run_gensa_for_seed <- function(seed) {
  opt_res <- GenSA(
    par = seed,
    fn = function(coord) market_share_objective(coord, zabka_medoids, POP.WAW, POW.waw, city_poly),
    lower = lower,
    upper = upper,
    control = list(max.call = 500, verbose = FALSE)
  )
  list(
    value = opt_res$value,
    par = opt_res$par
  )
}

# Parallel version: returns a list of results
cat("Running GenSA in parallel on", nrow(seed_coords), "seeds...\n")
results <- future_lapply(1:nrow(seed_coords), function(k) {
  seed <- seed_coords[k, ]
  res <- run_gensa_for_seed(seed)
  cat(sprintf("Seed %d: HHI %.8f\n", k, res$value))
  res
})

# Extract the best result
best_idx <- which.min(sapply(results, function(r) r$value))
best_hhi <- results[[best_idx]]$value
best_coord <- results[[best_idx]]$par
pt <- st_sfc(st_point(best_coord), crs = 2180)
best_result <- evaluate_candidate_impact(st_sf(geometry = pt), zabka_medoids, POP.WAW, POW.waw)
best_candidate <- st_sf(geometry = pt)

cat(sprintf("Best candidate found via optimization: HHI = %.5f\n", best_hhi))
cat("Coordinates (WGS84):\n")
print(st_coordinates(st_transform(best_candidate, 4326)))



# ---- 8. VISUALIZE CANDIDATES (BEFORE/AFTER) ----

plot_market_share_map_side_by_side <- function(vor_before, vor_after, school_pts, pow_sf, new_school_pt, river, railway) {
  # Calculate market shares & HHI
  vor_before$market_share <- vor_before$pop_tile / sum(vor_before$pop_tile, na.rm = TRUE)
  vor_after$market_share  <- vor_after$pop_tile  / sum(vor_after$pop_tile,  na.rm = TRUE)
  hhi_before <- sum(vor_before$spatially_weighted_share^2, na.rm = TRUE)
  hhi_after  <- sum(vor_after$spatially_weighted_share^2,  na.rm = TRUE)
  
  # Color scale
  all_shares <- c(vor_before$spatially_weighted_share, vor_after$spatially_weighted_share)
  n_breaks <- 6
  share_breaks <- seq(min(all_shares, na.rm = TRUE), max(all_shares, na.rm = TRUE), length.out = n_breaks + 1)
  share_colors <- colorRampPalette(c("yellow", "purple"))(n_breaks)
  vor_before$color <- cut(vor_before$spatially_weighted_share, breaks = share_breaks, labels = share_colors, include.lowest = TRUE)
  vor_after$color  <- cut(vor_after$spatially_weighted_share,  breaks = share_breaks, labels = share_colors, include.lowest = TRUE)
  par(mfrow = c(1, 2), mar = c(1, 1, 4, 1))
  # BEFORE
  plot(st_geometry(pow_sf), col = "white", border = "gray30",
       main = paste0("Before (current)\nHHI: ", round(hhi_before, 6)))
  plot(st_geometry(st_transform(vor_before, 4326)), col = as.character(vor_before$color), border = "black", lwd = 1.1, add = TRUE)
  plot(st_geometry(st_transform(st_as_sf(school_pts), 4326)), col = "black", bg = "green", pch = 21, cex = 1.3, add = TRUE)
  plot(st_geometry(river$osm_lines), col = "steelblue1", lwd = 3, add = TRUE)
  plot(st_geometry(railway$osm_lines), col = "gray", lwd = 1, add = TRUE)
  # AFTER
  plot(st_geometry(pow_sf), col = "white", border = "gray30",
       main = paste0("After (with new Zabka)\nHHI: ", round(hhi_after, 6)))
  plot(st_geometry(st_transform(vor_after, 4326)), col = as.character(vor_after$color), border = "black", lwd = 1.1, add = TRUE)
  plot(st_geometry(st_transform(st_as_sf(school_pts), 4326)), col = "black", bg = "green", pch = 21, cex = 1.3, add = TRUE)
  plot(st_geometry(st_transform(new_school_pt, 4326)), col = "red", bg = "white", pch = 8, cex = 1.0, lwd = 3, add = TRUE)
  plot(st_geometry(river$osm_lines), col = "steelblue1", lwd = 3, add = TRUE)
  plot(st_geometry(railway$osm_lines), col = "gray", lwd = 1, add = TRUE)
  par(xpd = NA)
  legend("bottomright",
         legend = paste0(round(share_breaks[-length(share_breaks)], 4), " to ", round(share_breaks[-1], 4)),
         fill = share_colors,
         title = "Market Share",
         bty = "n", cex = 0.8)
  par(mfrow = c(1, 1))
}

# ---- CALL THE FUNCTION ----

VOR_2180$pop_tile <- VOR_2180$pop_tile.x

sw_before <- compute_spatially_weighted_hhi(VOR_2180)
VOR_2180 <- sw_before$vor_sf

best_result$voronoi <- compute_spatially_weighted_hhi(best_result$voronoi)$vor_sf

plot_market_share_map_side_by_side(
  vor_before = VOR_2180, 
  vor_after  = best_result$voronoi,  
  school_pts = zabka_medoids, 
  pow_sf     = POW.waw, 
  new_school_pt = best_candidate, 
  river = river, 
  railway = railway
)

