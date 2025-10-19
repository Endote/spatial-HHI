# ========================== SPATIALLY WEIGHTED HHI: LIDL vs BIEDRONKA (WARSZAWA) ==========================
# RUN-READY SINGLE FILE — STATIC BREAKS (same as ETAP 2)

# ---- Packages ----
req <- c("sf","dplyr","osmdata","RColorBrewer","spatstat.geom","units","igraph")
invisible(lapply(req, function(p) if (!requireNamespace(p, quietly=TRUE)) install.packages(p)))
lapply(req, library, character.only = TRUE)
options(scipen = 999)

# ---- CONFIG ----
use_local_boundary <- FALSE     # TRUE jeśli masz powiat/dzielnice z pliku
use_local_popgrid  <- TRUE      # TRUE jeśli masz gotowe POP.WAW z kolumną RES (zalecane)
alpha_NSCI <- 0.6               # waga między-sieciowa w NSCI(α)
random_grid_total_pop <- 1e6    # fallback: łączna ludność na jednorodnej siatce

# === STATIC BREAKS (identyczne jak w ETAP 2) ===
S_breaks_global <- c(0.000, 0.004, 0.008, 0.012, 0.016, 0.030)

# ---- ŚCIEŻKI OPCJONALNYCH WARSTW LOKALNYCH ----
path_powiaty   <- "data/powiaty/powiaty.shp"
path_dzielnice <- "data/dzielnice_Warszawy/dzielnice_Warszawy.shp"
path_pop_waw   <- "data/POP.WAW/POP.WAW.shp"

# ---- Helper: scalanie klastrów (medoid) ----
merge_points_medoid <- function(points_sf, threshold_m = 200, metric_crs = 2180, output_crs = 4326) {
  if (nrow(points_sf) == 0) return(points_sf)
  pts_m <- st_transform(points_sf, metric_crs)
  dm <- st_distance(pts_m)
  thr <- set_units(threshold_m, "m")
  adj <- (dm < thr) & (dm > set_units(0,"m"))
  g <- graph_from_adjacency_matrix(as.matrix(adj), mode="undirected", diag=FALSE)
  comp <- components(g)
  if (length(comp$membership) == 0) {
    out <- pts_m
  } else {
    pts_m$cluster <- comp$membership
    medoid_idx <- c()
    for (cl in sort(unique(pts_m$cluster))) {
      idx <- which(pts_m$cluster == cl)
      if (length(idx) == 1) {
        medoid_idx <- c(medoid_idx, idx)
      } else {
        sub_dm <- as.matrix(dm[idx, idx])
        medoid_idx <- c(medoid_idx, idx[which.min(rowSums(sub_dm))])
      }
    }
    out <- pts_m[unique(medoid_idx), ]
  }
  st_transform(out, output_crs)
}

# ---- Boundary (OSM fallback) ----
message("Fetching Warsaw boundary…")
if (use_local_boundary && file.exists(path_powiaty)) {
  POW <- st_read(path_powiaty, quiet = TRUE) |> st_transform(4326)
  names(POW) <- iconv(names(POW), from="latin1", to="UTF-8")
  nmcol <- grep("jpt_nazwa", names(POW), value = TRUE)[1]
  POW.waw <- POW[POW[[nmcol]] %in% c("powiat Warszawa","Warszawa"),]
  if (nrow(POW.waw)==0) stop("Nie znaleziono 'powiat Warszawa' w warstwie powiatów.")
  city_boundary_4326 <- st_union(POW.waw)
} else {
  bb <- getbb("Warsaw, Poland")
  adm <- opq(bbox = bb) |>
    add_osm_feature(key="admin_level", value=c("6","7","8")) |>
    add_osm_feature(key="name", value="Warszawa") |>
    osmdata_sf()
  cand <- list(adm$osm_multipolygons, adm$osm_polygons)
  cand <- cand[!sapply(cand, is.null)]
  if (length(cand)==0) stop("Brak granicy Warszawy z OSM.")
  geoms <- do.call(rbind, lapply(cand, function(x) st_cast(st_make_valid(x), "MULTIPOLYGON")))
  geoms <- geoms[st_is_valid(geoms),]
  geoms$area <- st_area(geoms)
  city_boundary_4326 <- geoms[which.max(geoms$area),] |> st_geometry() |> st_union()
}
city_boundary_2180 <- st_transform(city_boundary_4326, 2180)

# Dzielnice (opcjonalne) – używamy tylko do clipowania punktów; fallback = boundary
WAW <- if (use_local_boundary && file.exists(path_dzielnice)) {
  st_read(path_dzielnice, quiet = TRUE) |> st_transform(4326)
} else {
  st_as_sf(st_make_valid(city_boundary_4326))
}

# ---- Kontekst (opcjonalny): rzeki i kolej ----
message("Fetching rivers & railways…")
warsaw_bbox <- st_bbox(city_boundary_4326)
river <- opq(bbox = warsaw_bbox) |> add_osm_feature(key="waterway", value="river") |> osmdata_sf()
railway <- opq(bbox = warsaw_bbox) |> add_osm_feature(key="railway", value="rail") |> osmdata_sf()

# ---- Punkty OSM: Biedronka i Lidl ----
message("Fetching Biedronka & Lidl points…")
q_b <- opq(bbox = warsaw_bbox) |> add_osm_feature(key="shop", value="supermarket") |> add_osm_feature(key="name", value="Biedronka")
q_l <- opq(bbox = warsaw_bbox) |> add_osm_feature(key="shop", value="supermarket") |> add_osm_feature(key="name", value="Lidl")
b_sf <- try(osmdata_sf(q_b)$osm_points, silent=TRUE)
l_sf <- try(osmdata_sf(q_l)$osm_points, silent=TRUE)
if (inherits(b_sf, "try-error") || inherits(l_sf, "try-error")) stop("Problem z pobraniem danych OSM.")

b_sf <- suppressWarnings(st_intersection(st_make_valid(b_sf), st_make_valid(WAW)))
l_sf <- suppressWarnings(st_intersection(st_make_valid(l_sf), st_make_valid(WAW)))
b_sf$supermarket <- "Biedronka"; l_sf$supermarket <- "Lidl"

biedronka_medoids <- merge_points_medoid(b_sf, threshold_m=200)
lidl_medoids      <- merge_points_medoid(l_sf, threshold_m=200)

supermarket_points <- dplyr::bind_rows(biedronka_medoids, lidl_medoids) |> st_transform(4326)
if (nrow(supermarket_points) == 0) stop("Brak punktów supermarketów po czyszczeniu.")

# ---- Siatka populacji (preferowana: POP.WAW). Fallback: jednorodna siatka ----
message("Preparing population grid…")
if (use_local_popgrid && file.exists(path_pop_waw)) {
  POP_WAW <- st_read(path_pop_waw, quiet = TRUE) |> st_transform(4326)
  if (!("RES" %in% names(POP_WAW))) stop("POP.WAW bez kolumny 'RES'.")
} else {
  city_2180 <- st_transform(city_boundary_4326, 2180)
  sz <- 800
  g <- st_make_grid(city_2180, cellsize = c(sz, sz), what="polygons", square=TRUE)
  POP_WAW <- st_sf(geometry = st_intersection(g, city_2180)) |> st_transform(4326)
  POP_WAW$RES <- random_grid_total_pop / nrow(POP_WAW)
}
POP_2180 <- st_transform(POP_WAW, 2180) |> mutate(cell_id = row_number(), cell_area = st_area(geometry))

# ---- Voronoi ALL (wszystkie punkty) + clip jako sf ----
message("Building Voronoi (ALL)…")
pts_2180 <- st_transform(supermarket_points, 2180)
env <- st_union(city_boundary_2180)
env_sf <- st_sf(geometry = env); st_crs(env_sf) <- 2180

vor_sfc <- st_voronoi(st_union(st_geometry(pts_2180)), envelope = env)
vor_sf  <- st_sf(geometry = st_collection_extract(vor_sfc, "POLYGON")); st_crs(vor_sf) <- 2180
VOR_2180 <- suppressWarnings(st_intersection(vor_sf, env_sf)) |> st_make_valid()
VOR_2180$tile_id  <- seq_len(nrow(VOR_2180))
VOR_2180$tile_area <- st_area(VOR_2180)

# przypis sieci do poligonów; fallback: najbliższy punkt
VOR_2180 <- st_join(VOR_2180, st_transform(supermarket_points, 2180) %>% select(supermarket), join = st_contains)
if (any(is.na(VOR_2180$supermarket))) {
  nearest_idx <- st_nearest_feature(st_centroid(VOR_2180), st_transform(supermarket_points, 2180))
  VOR_2180$supermarket[is.na(VOR_2180$supermarket)] <- st_transform(supermarket_points, 2180)$supermarket[nearest_idx[is.na(VOR_2180$supermarket)]]
}

# ---- POP ∩ VOR (ALL) – frakcyjne przypisanie populacji ----
message("Intersecting population with Voronoi (ALL)…")
pop_vor <- suppressWarnings(st_intersection(POP_2180, VOR_2180)) |>
  mutate(sub_area = st_area(geometry),
         pop_share = RES * (sub_area / pmax(cell_area, set_units(1,"m^2"))))

pop_per_tile <- pop_vor |>
  st_drop_geometry() |>
  group_by(tile_id) |>
  summarise(pop_tile = sum(pop_share, na.rm=TRUE), .groups="drop")

VOR_2180 <- VOR_2180 |>
  left_join(pop_per_tile, by="tile_id") |>
  mutate(pop_tile = ifelse(is.na(pop_tile), 0, as.numeric(pop_tile)),
         cell_area_km2 = as.numeric(tile_area) / 1e6,
         pop_density = ifelse(cell_area_km2>0, pop_tile / cell_area_km2, 0))

# ---- Funkcje liczące S_i/HHI ----
compute_Si_HHI <- function(vor_subset) {
  Pi <- as.numeric(vor_subset$pop_tile)
  Ai <- as.numeric(st_area(vor_subset))
  Pi[is.na(Pi)] <- 0
  Ai[is.na(Ai) | Ai==0] <- 1e-9
  pi  <- Pi / Ai
  wi  <- if (sum(pi, na.rm=TRUE)>0) pi / sum(pi, na.rm=TRUE) else rep(0, length(pi))
  Ptot <- sum(Pi, na.rm=TRUE)
  ms  <- if (Ptot > 0) Pi / Ptot else rep(0, length(Pi))
  S_i <- ms * wi
  S_i <- if (sum(S_i, na.rm=TRUE) > 0) S_i / sum(S_i, na.rm=TRUE) else S_i
  list(S_i=S_i, HHI_ALL=sum(S_i^2, na.rm=TRUE))
}

# (NOWOŚĆ) – pełny pipeline INTRA dla wybranej sieci: własny Voronoi + POP ∩ Vor + S_i
compute_chain_voronoi <- function(points_chain_4326, boundary_2180, POP_2180) {
  if (is.null(points_chain_4326) || nrow(points_chain_4326) == 0) {
    return(list(vor = st_sf(geometry = st_sfc(), crs = 2180), HHI = NA_real_))
  }
  pts <- st_transform(points_chain_4326, 2180)
  env <- st_union(boundary_2180)
  env_sf <- st_sf(geometry = env); st_crs(env_sf) <- st_crs(boundary_2180)
  
  vor_sfc <- st_voronoi(st_union(st_geometry(pts)), envelope = env)
  vor_sf  <- st_sf(geometry = st_collection_extract(vor_sfc, "POLYGON")); st_crs(vor_sf) <- st_crs(boundary_2180)
  vor_sf  <- suppressWarnings(st_intersection(vor_sf, env_sf)) |> st_make_valid()
  vor_sf$tile_id <- seq_len(nrow(vor_sf))
  
  pop_vor <- suppressWarnings(st_intersection(POP_2180, vor_sf)) |>
    mutate(sub_area = st_area(geometry),
           pop_share = RES * (sub_area / pmax(cell_area, set_units(1, "m^2"))))
  
  pop_per_tile <- pop_vor |>
    st_drop_geometry() |>
    group_by(tile_id) |>
    summarise(pop_tile = sum(pop_share, na.rm = TRUE), .groups = "drop")
  
  vor_sf <- vor_sf |>
    left_join(pop_per_tile, by = "tile_id") |>
    mutate(pop_tile = ifelse(is.na(pop_tile), 0, as.numeric(pop_tile)))
  
  Pi <- as.numeric(vor_sf$pop_tile)
  Ai <- as.numeric(st_area(vor_sf)); Ai[Ai<=0 | is.na(Ai)] <- 1e-9
  pi <- Pi / Ai
  wi <- if (sum(pi, na.rm=TRUE)>0) pi / sum(pi, na.rm=TRUE) else rep(0, length(pi))
  Ptot <- sum(Pi, na.rm=TRUE)
  ms <- if (Ptot > 0) Pi / Ptot else rep(0, length(Pi))
  S_i <- ms * wi
  if (sum(S_i, na.rm=TRUE) > 0) S_i <- S_i / sum(S_i, na.rm=TRUE)
  
  vor_sf$S_i_intra <- S_i
  HHI_intra <- sum(S_i^2, na.rm = TRUE)
  list(vor = vor_sf, HHI = HHI_intra)
}

# ---- ALL: S_i i HHI_ALL + HHI_INTER ----
res_all <- compute_Si_HHI(VOR_2180)
VOR_2180$S_i_all <- res_all$S_i
HHI_ALL <- res_all$HHI_ALL

chain_shares <- VOR_2180 |>
  st_drop_geometry() |>
  group_by(supermarket) |>
  summarise(S_c = sum(S_i_all, na.rm=TRUE), .groups="drop")
HHI_INTER <- sum(chain_shares$S_c^2, na.rm=TRUE)

# ---- INTRA: oddzielne Voronoi i POP ∩ Vor dla każdej sieci ----
b_pts_4326 <- supermarket_points %>% filter(supermarket == "Biedronka")
l_pts_4326 <- supermarket_points %>% filter(supermarket == "Lidl")
res_B <- compute_chain_voronoi(b_pts_4326, city_boundary_2180, POP_2180)
res_L <- compute_chain_voronoi(l_pts_4326, city_boundary_2180, POP_2180)

# ---- NSCI(α) ----
NSCI_alpha <- alpha_NSCI * HHI_INTER + (1 - alpha_NSCI) * HHI_ALL
message(sprintf("HHI_ALL = %.3f, HHI_INTER = %.3f, NSCI(α=%.2f) = %.3f", HHI_ALL, HHI_INTER, alpha_NSCI, NSCI_alpha))


# ---- Binner z STAŁYMI progami (jak w ETAP 2) ----
# używa findInterval, więc wartości < min trafiają do 1 binu, > max do ostatniego binu
make_binned_palette <- function(vals, breaks) {
  v_all <- as.numeric(vals)
  br <- sort(unique(as.numeric(breaks)))
  if (length(v_all) == 0 || length(br) < 2) {
    return(list(cols = rep("gray90", length(v_all)),
                legend_labels = "no data",
                legend_cols = "gray90",
                breaks = br))
  }
  nb <- length(br) - 1
  pal <- colorRampPalette(c("yellow","purple"))(nb)
  # indeksy binów 1..nb (clampowane do 1..nb)
  idx <- findInterval(v_all, vec = br, rightmost.closed = TRUE, all.inside = TRUE)
  idx[idx < 1]  <- 1
  idx[idx > nb] <- nb
  cols_vec <- pal[idx]
  lab <- paste0(format(br[-length(br)], digits = 3, nsmall = 3),
                " – ",
                format(br[-1],          digits = 3, nsmall = 3))
  list(cols = cols_vec, legend_labels = lab, legend_cols = pal, breaks = br)
}

# ---- Rysunki: 1×3 (ALL | INTRA Biedronka | INTRA Lidl) ----
par(mfrow=c(1,3), mar=c(1.5,2,3.6,1.5), oma=c(0,0,2,0))

# 1) ALL — STAŁE progi S_breaks_global
bins_all <- make_binned_palette(VOR_2180$S_i_all, breaks = S_breaks_global)
plot(st_transform(city_boundary_4326, 4326), col="white", border="gray30",
     main = sprintf("ALL: HHI_ALL = %.3f,  HHI_INTER = %.3f", HHI_ALL, HHI_INTER), reset=FALSE)
plot(st_transform(VOR_2180, 4326) |> st_geometry(), col = bins_all$cols, border="black", add=TRUE)
pts_4326 <- st_transform(supermarket_points, 4326)
plot(st_geometry(pts_4326), pch=21, cex=1.5, col="black",
     bg = ifelse(pts_4326$supermarket=="Biedronka","red","blue"), add=TRUE)
legend("bottomleft",
       legend = bins_all$legend_labels,
       fill   = bins_all$legend_cols,
       border = "gray30",
       title  = "S_i (ALL)",
       bty    = "n", cex = 0.8)
if (!is.null(river$osm_lines))   plot(st_geometry(river$osm_lines),   col="lightblue", lwd=2, add=TRUE)
if (!is.null(railway$osm_lines)) plot(st_geometry(railway$osm_lines), col="gray70",   lwd=1, add=TRUE)
# liczba obserwacji (nie-NA)
n_all <- sum(!is.na(VOR_2180$S_i_all))
mtext(sprintf("Number of observations: %d", n_all), side=3, adj=1, line=-1.5, cex=0.8)

# 2) INTRA Biedronka — STAŁE progi S_breaks_global
plot(st_transform(city_boundary_4326, 4326), col="white", border="gray30",
     main = sprintf("INTRA Biedronka: HHI_B = %.3f", res_B$HHI), reset=FALSE)
if (nrow(res_B$vor) > 0) {
  bins_B <- make_binned_palette(res_B$vor$S_i_intra, breaks = S_breaks_global)
  plot(st_geometry(st_transform(res_B$vor, 4326)), col=bins_B$cols, border="black", add=TRUE)
  b_pts <- st_transform(b_pts_4326, 4326)
  plot(st_geometry(b_pts), pch=21, cex=1.5, col="black", bg="red", add=TRUE)
  legend("bottomleft",
         legend = bins_B$legend_labels,
         fill   = bins_B$legend_cols,
         border = "gray30",
         title  = "S_i (INTRA B)",
         bty    = "n", cex = 0.8)
  n_B <- sum(!is.na(res_B$vor$S_i_intra))
  mtext(sprintf("Number of observations: %d", n_B), side=3, adj=1, line=-1.5, cex=0.8)
} else {
  mtext("Brak sklepów Biedronka w obszarze", side=3, line=-1)
}
if (!is.null(river$osm_lines))   plot(st_geometry(river$osm_lines),   col="lightblue", lwd=2, add=TRUE)
if (!is.null(railway$osm_lines)) plot(st_geometry(railway$osm_lines), col="gray70",   lwd=1, add=TRUE)

# 3) INTRA Lidl — STAŁE progi S_breaks_global
plot(st_transform(city_boundary_4326, 4326), col="white", border="gray30",
     main = sprintf("INTRA Lidl: HHI_L = %.3f", res_L$HHI), reset=FALSE)
if (nrow(res_L$vor) > 0) {
  bins_L <- make_binned_palette(res_L$vor$S_i_intra, breaks = S_breaks_global)
  plot(st_geometry(st_transform(res_L$vor, 4326)), col=bins_L$cols, border="black", add=TRUE)
  l_pts <- st_transform(l_pts_4326, 4326)
  plot(st_geometry(l_pts), pch=21, cex=1.5, col="black", bg="blue", add=TRUE)
  legend("bottomleft",
         legend = bins_L$legend_labels,
         fill   = bins_L$legend_cols,
         border = "gray30",
         title  = "S_i (INTRA L)",
         bty    = "n", cex = 0.8)
  n_L <- sum(!is.na(res_L$vor$S_i_intra))
  mtext(sprintf("Number of observations: %d", n_L), side=3, adj=1, line=-1.5, cex=0.8)
} else {
  mtext("Brak sklepów Lidl w obszarze", side=3, line=-1)
}
if (!is.null(river$osm_lines))   plot(st_geometry(river$osm_lines),   col="lightblue", lwd=2, add=TRUE)
if (!is.null(railway$osm_lines)) plot(st_geometry(railway$osm_lines), col="gray70",   lwd=1, add=TRUE)

mtext(sprintf("Lidl & Biedronka — Spatially Weighted HHI | NSCI(α=%.2f)=%.3f",
              alpha_NSCI, NSCI_alpha),
      side=3, outer=TRUE, line=0.5, cex=1.1, font=2)

# ========================== END ==========================
