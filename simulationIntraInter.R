# --- Libraries ---
libs <- c("sf", "dplyr")
invisible(lapply(libs, function(pkg) if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)))
lapply(libs, library, character.only = TRUE)

# ----------------- PARAMS -----------------
set.seed(123)                       # powtarzalność
n_per_network <- c(A=4, B=2, C=3)   # <<<< 4 / 2 / 3
n_pop_side <- 32
total_pop <- 5000
crs_geo  <- 4326
crs_proj <- 3857

# ----------------- BOUNDARY & POP -----------------
boundary <- st_polygon(list(rbind(c(0,0), c(1,0), c(1,1), c(0,1), c(0,0)))) |>
  st_sfc(crs = crs_geo)
boundary_proj <- st_transform(boundary, crs_proj)

pop_tiles <- st_make_grid(boundary, n = c(n_pop_side, n_pop_side)) |> st_sf()
pop_tiles_proj <- st_transform(pop_tiles, crs_proj)
n_tiles <- nrow(pop_tiles_proj)
pop_tiles_proj$pop <- rep(total_pop / n_tiles, n_tiles)
pop_tiles_proj$tile_id <- 1:n_tiles

# ----------------- LOSOWE PLACÓWKI 3 SIECI -----------------
rand_pts <- function(n) data.frame(x = runif(n, 0, 1), y = runif(n, 0, 1))
pts_df <- dplyr::bind_rows(
  cbind(rand_pts(n_per_network["A"]), network="A"),
  cbind(rand_pts(n_per_network["B"]), network="B"),
  cbind(rand_pts(n_per_network["C"]), network="C")
)
pts_sf  <- st_as_sf(pts_df, coords=c("x","y"), crs = crs_geo) |> st_transform(crs_proj)
pts_sf$store_id <- 1:nrow(pts_sf)

# ----------------- FUNKCJE -----------------
compute_shares <- function(points_sf, pop_tiles_proj, boundary_proj){
  vor <- st_voronoi(st_union(points_sf), envelope = boundary_proj)
  vor_polys <- st_collection_extract(vor, "POLYGON")
  vor_sf <- st_sf(geometry = vor_polys)
  vor_sf <- st_intersection(vor_sf, boundary_proj)
  vor_sf$store_id <- points_sf$store_id
  
  pop_vor <- suppressWarnings(st_intersection(pop_tiles_proj, vor_sf))
  pop_vor$sub_area <- as.numeric(st_area(pop_vor))
  tile_areas <- as.numeric(st_area(pop_tiles_proj))
  names(tile_areas) <- pop_tiles_proj$tile_id
  pop_vor$tile_area <- tile_areas[as.character(pop_vor$tile_id)]
  pop_vor$pop_share <- pop_vor$pop * (pop_vor$sub_area / pmax(pop_vor$tile_area, 1e-9))
  
  pop_sum <- pop_vor |>
    st_drop_geometry() |>
    group_by(store_id) |>
    summarize(pop_tile = sum(pop_share, na.rm = TRUE), .groups='drop')
  vor_sf <- left_join(vor_sf, pop_sum, by = "store_id")
  vor_sf$pop_tile[is.na(vor_sf$pop_tile)] <- 0
  
  Pi <- vor_sf$pop_tile
  Ai <- as.numeric(st_area(vor_sf))
  pi <- Pi / pmax(Ai, 1e-9)
  wi <- pi / sum(pi, na.rm=TRUE)
  Ptot <- sum(Pi, na.rm=TRUE)
  ms <- if (Ptot > 0) Pi / Ptot else rep(0, length(Pi))
  S_i <- ms * wi
  S_i <- S_i / sum(S_i, na.rm=TRUE)
  
  vor_sf$S_i <- S_i
  list(vor = vor_sf, S_i = S_i, HHI_all = sum(S_i^2, na.rm=TRUE))
}

palette_from_values <- function(vals){
  vals_r <- round(vals, 4)
  u <- sort(unique(vals_r))
  pal <- colorRampPalette(c("yellow","purple"))(length(u))
  col <- pal[match(vals_r, u)]
  list(cols = col, legend_vals = u, legend_cols = pal)
}

# ----------------- OBLICZENIA GLOBALNE (ALL) -----------------
res_all <- compute_shares(pts_sf, pop_tiles_proj, boundary_proj)
vor_all <- res_all$vor
HHI_ALL <- res_all$HHI_all

network_shares <- vor_all |>
  st_drop_geometry() |>
  left_join(st_drop_geometry(pts_sf[,c("store_id","network")]), by="store_id") |>
  group_by(network) |>
  summarize(S_c = sum(S_i, na.rm=TRUE), .groups='drop')
HHI_INTER <- sum(network_shares$S_c^2)

# ----------------- HHI wewnątrz sieci -----------------
compute_intra_for_network <- function(net_label){
  pts_net <- pts_sf[pts_sf$network==net_label,]
  pts_net$store_id <- 1:nrow(pts_net)
  res <- compute_shares(pts_net, pop_tiles_proj, boundary_proj)
  list(vor=res$vor, HHI_INTRA=res$HHI_all, pts=pts_net)
}
res_A <- compute_intra_for_network("A")
res_B <- compute_intra_for_network("B")
res_C <- compute_intra_for_network("C")

# ----------------- RYSUNKI -----------------
par(mfrow=c(2,2), mar=c(1.5,2,3.6,1.5), oma=c(0,0,2,0))

cols_all <- palette_from_values(vor_all$S_i)
plot(st_geometry(vor_all), col = cols_all$cols, border="gray30")
points(st_coordinates(pts_sf),
       pch = c(8,24,21)[as.integer(factor(pts_sf$network, levels=c("A","B","C")))],
       cex = 1.3, bg = "white")
legend("topright", inset = c(0.02,0.02),
       legend = sprintf("%.3f", cols_all$legend_vals),
       fill = cols_all$legend_cols, border="gray30", title = "S_i", bty="n", cex=0.9)
mtext(sprintf("ALL: HHI_ALL = %.3f,  HHI_INTER = %.3f", HHI_ALL, HHI_INTER),
      side=3, line=0.5, cex=0.95, font=2)

cols_A <- palette_from_values(res_A$vor$S_i)
plot(st_geometry(res_A$vor), col = cols_A$cols, border="gray30")
points(st_coordinates(res_A$pts), pch=8, cex=1.3, bg="white")
mtext(sprintf("INTRA A (4 stores): HHI_A = %.3f", res_A$HHI_INTRA),
      side=3, line=0.5, cex=0.95, font=2)

cols_B <- palette_from_values(res_B$vor$S_i)
plot(st_geometry(res_B$vor), col = cols_B$cols, border="gray30")
points(st_coordinates(res_B$pts), pch=24, cex=1.3, bg="white")
mtext(sprintf("INTRA B (2 stores): HHI_B = %.3f", res_B$HHI_INTRA),
      side=3, line=0.5, cex=0.95, font=2)

cols_C <- palette_from_values(res_C$vor$S_i)
plot(st_geometry(res_C$vor), col = cols_C$cols, border="gray30")
points(st_coordinates(res_C$pts), pch=21, cex=1.3, bg="white")
mtext(sprintf("INTRA C (3 stores): HHI_C = %.3f", res_C$HHI_INTRA),
      side=3, line=0.5, cex=0.95, font=2)

mtext("Three networks: ALL vs INTRA — counts A=4, B=2, C=3",
      side=3, outer=TRUE, line=0.5, cex=1.1, font=2)

# --- (opcjonalnie) indeks mieszany NSCI(α) ---
alpha <- 0.6
NSCI_alpha <- alpha*HHI_INTER + (1-alpha)*HHI_ALL
cat(sprintf("\nHHI_ALL = %.4f,  HHI_INTER = %.4f,  NSCI(α=%.2f) = %.4f\n",
            HHI_ALL, HHI_INTER, alpha, NSCI_alpha))
