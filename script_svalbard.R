#### MAP ####
# Map of sampling sites with bathymetry and GPS tracks - Svalbard
# This section produces a high-resolution map showing sampling sites overlaid
# on bathymetry and coastline data for the Svalbard region.

# Libraries
library(sf)                 # Handling spatial vector data
library(ggplot2)            # Plotting framework
library(rnaturalearth)      # Natural Earth map data
library(marmap)             # Bathymetry data from NOAA
library(viridis)            # Colorblind-friendly color scales
library(readr)              # Reading tabular data
library(ggrepel)            # Non-overlapping text labels
library(dplyr)              # Data manipulation
library(rnaturalearthdata)  # Supporting data for rnaturalearth

# Load sampling site data (columns: Site, Category, Latitude, Longitude)
data <- read_csv2("svalbard_sites.csv")
data <- data %>% rename(Category = Catégorie)
# Convert to an sf object using WGS84 coordinates
data_sf <- st_as_sf(data, coords = c("Longitude", "Latitude"), crs = 4326)

# Bathymetry limits (adjust if needed)
# Download NOAA bathymetry for the study area
bathy <- getNOAA.bathy(lon1 = 10, lon2 = 20, lat1 = 78, lat2 = 80.5, resolution = 0.05)

# Convert bathymetry to data frame for ggplot
bathy_df <- marmap::fortify.bathy(bathy)

# Use only negative depths (seafloor)
bathy_df_neg <- bathy_df %>% filter(z < 0)

# Svalbard bounding box
bbox_svalbard <- st_as_sfc(st_bbox(c(xmin = 10, xmax = 20, ymin = 78, ymax = 80.5), crs = 4326))

# Land shapefiles
# Load country polygons and crop to the study area
land <- ne_countries(scale = "medium", returnclass = "sf")
land_crop <- st_crop(land, bbox_svalbard)

# High-resolution coastline
# Load and crop coastline geometry
coastline <- ne_download(scale = 10, type = "coastline", category = "physical", returnclass = "sf")
coastline_crop <- st_crop(coastline, bbox_svalbard)

# Plot
p <- ggplot() +
  
  # Bathymetry raster layer
  geom_raster(data = bathy_df_neg, aes(x = x, y = y, fill = z)) +
  scale_fill_viridis(
    option = "mako",
    direction = 1,
    name = "Depth (m)",
    limits = c(min(bathy_df_neg$z), 0),
    oob = scales::squish
  ) +
  
  # Depth contour lines
  geom_contour(
    data = bathy_df_neg,
    aes(x = x, y = y, z = z),
    breaks = c(-200, -1000, -2000),
    color = "grey",
    size = 0.3
  ) +
  
  # Sampling sites
  geom_sf(data = data_sf, aes(shape = Category), size = 2) +
  
  # Site labels
  # Repelled text labels to avoid overlap
  geom_text_repel(
    data = data_sf,
    aes(
      x = st_coordinates(data_sf)[,1],
      y = st_coordinates(data_sf)[,2],
      label = Site
    ),
    color = "black",
    size = 3.5,
    box.padding = 0.5,
    segment.color = "grey40",
    segment.size = 0.3,
    bg.color = "white",
    bg.r = 0.15
  ) +
  
  # Shapes for site categories
  scale_shape_manual(
    values = c(15, 16, 17, 15, 18),
    name = "Site category"
  ) +
  
  # Map extent
  coord_sf(xlim = c(10, 16), ylim = c(78, 80), expand = FALSE) +
  
  # Axis labels and theme
  labs(
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10)
  )

# Print the plot
p

# Export the map as a high-resolution PNG
ggsave(
  filename = "svalbard_map_highres.png",
  plot = p,
  width = 8,        # width in inches
  height = 6,       # height in inches
  dpi = 600
)


#### DATA ####
# Data aggregation and transformation for eDNA presence/absence analysis

library(dplyr)
library(tidyr)
library(stringr)

# Read input files
data <- read.csv2("data.csv", stringsAsFactors = FALSE, check.names = FALSE)
meta <- read.csv2("metadata.csv", stringsAsFactors = FALSE, check.names = FALSE)

# Columns starting with MUSA (sample identifiers)
musa_cols <- grep("^MUSA", names(data), value = TRUE)

# Convert read counts to binary presence/absence (1 = present, 0 = absent)
data_bin <- data %>%
  mutate(across(all_of(musa_cols), ~ as.integer(. > 0)))

# Extract everything before the first dash in the client ID
meta <- meta %>%
  mutate(ID_prefix = sub("-.*", "", `ID client`))

# Create a mapping between MUSA IDs and site prefixes
musa_to_prefix <- meta %>%
  select(`ID Argaly`, ID_prefix)

# Convert to long format, aggregate by prefix, and return to wide format
final_df <- data_bin %>%
  pivot_longer(cols = all_of(musa_cols), names_to = "MUSA", values_to = "val") %>%
  left_join(musa_to_prefix, by = c("MUSA" = "ID Argaly")) %>%
  group_by(across(-c(MUSA, val, ID_prefix)), ID_prefix) %>%
  summarise(val = as.integer(any(val == 1)), .groups = "drop") %>%
  pivot_wider(names_from = ID_prefix, values_from = val, values_fill = 0) %>%
  rename(List_E = `List E (a: possible taxon absent from Genbank, b: similar taxon but not distributed)`)

# Renaming vector: old site code = new standardized site name
prefix_renames <- c(
  "SFN" = "DEEP_SFN",
  "POO" = "WAL_POO",
  "SAR" = "SHA_SAR",
  "GUL" = "WAL_GUL",
  "SCH" = "SHA_SCH",
  "FUG" = "GLA_FUG",
  "SME1" = "WAL_SME",
  "SME2" = "GLA_SME", # Warning: duplicates
  "LIL" = "GLA_LIL",
  "NYA" = "HAR_NYA",
  "KON1" = "GLA_KON",  # Warning: duplicates
  "KON2" = "DEEP_KON",
  "DAH" = "GLA_DAH",
  "YME" = "GLA_YME",
  "BAR" = "HAR_BAR",
  "LON" = "HAR_LON"
)

# Get names of columns 19 to 34 (site columns)
old_names <- names(final_df)[19:34]

# Keep only columns that exist in the data frame
rename_map <- prefix_renames[names(prefix_renames) %in% old_names]

# Reverse mapping: new_name = old_name
rename_map <- setNames(names(rename_map), rename_map)

# Apply renaming and filter to target vertebrate classes
final_df <- final_df %>%
  rename(!!!rename_map) %>%
  filter(Class %in% c("Mammalia", "Aves", "Teleostei", "Elasmobranchii"))

# Export aggregated data
write.csv(final_df, "data_aggrege.csv", row.names = FALSE, fileEncoding = "UTF-8")


#### VENN ####
# Venn diagram showing taxonomic overlap among site categories

# Install the package if needed
# install.packages("ggvenn")

library(ggvenn)
library(dplyr)

# Columns corresponding to individual sampling sites
site_cols <- c("DEEP_SFN", "DEEP_KON", "SHA_SAR", "SHA_SCH", "GLA_FUG", 
               "GLA_SME", "GLA_LIL", "GLA_KON", "GLA_DAH", "GLA_YME", 
               "HAR_NYA", "HAR_BAR", "HAR_LON", "WAL_POO", "WAL_GUL", "WAL_SME")

# Create a list of taxa detected at each site
taxon_lists <- lapply(final_df[site_cols], function(x) final_df$Taxonomie[x == 1])

# Group taxa by environmental category
taxon_lists_grouped <- list(
  DEEP = unique(c(taxon_lists$DEEP_SFN, taxon_lists$DEEP_KON)),
  SHA = unique(c(taxon_lists$SHA_SAR, taxon_lists$SHA_SCH)),
  GLA = unique(c(taxon_lists$GLA_FUG, taxon_lists$GLA_SME, taxon_lists$GLA_LIL, taxon_lists$GLA_KON, taxon_lists$GLA_DAH, taxon_lists$GLA_YME)),
  HAR = unique(c(taxon_lists$HAR_NYA, taxon_lists$HAR_BAR, taxon_lists$HAR_LON)),
  WAL = unique(c(taxon_lists$WAL_POO, taxon_lists$WAL_GUL, taxon_lists$WAL_SME))
)

# Plot the Venn diagram
ggvenn(
  taxon_lists_grouped,
  fill_color = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF", "#FF61CC"),
  set_name_size = 5,   # size of category labels
  text_size = 5        # size of numbers in each intersection
)


#### TAXON LISTS, UNIQUENESS, INTERSECTIONS, AND BARPLOTS ####
# This section derives taxon lists per site and per environmental category,
# identifies taxa unique to each category, computes shared taxa among all
# possible category combinations, and visualizes taxonomic richness by class.

# List of taxa by site (same logic as above)
# Each element contains the taxa detected (presence = 1) at a given site
site_cols <- c("DEEP_SFN", "DEEP_KON", "SHA_SAR", "SHA_SCH", "GLA_FUG", 
               "GLA_SME", "GLA_LIL", "GLA_KON", "GLA_DAH", "GLA_YME", 
               "HAR_NYA", "HAR_BAR", "HAR_LON", "WAL_POO", "WAL_GUL", "WAL_SME")

taxon_lists <- lapply(final_df[site_cols], function(x) final_df$Taxonomie[x == 1])

# Group taxa by environmental category
taxon_lists_grouped <- list(
  DEEP = unique(c(taxon_lists$DEEP_SFN, taxon_lists$DEEP_KON)),
  SHA  = unique(c(taxon_lists$SHA_SAR,  taxon_lists$SHA_SCH)),
  GLA  = unique(c(taxon_lists$GLA_FUG, taxon_lists$GLA_SME, taxon_lists$GLA_LIL,
                  taxon_lists$GLA_KON, taxon_lists$GLA_DAH, taxon_lists$GLA_YME)),
  HAR  = unique(c(taxon_lists$HAR_NYA, taxon_lists$HAR_BAR, taxon_lists$HAR_LON)),
  WAL  = unique(c(taxon_lists$WAL_POO, taxon_lists$WAL_GUL, taxon_lists$WAL_SME))
)

# Alternative explicit formulation (same result, more readable)
# Lists of taxa per category extracted directly from the data frame
taxon_lists_grouped <- list(
  DEEP = unique(c(final_df$Taxonomie[final_df$DEEP_SFN == 1],
                  final_df$Taxonomie[final_df$DEEP_KON == 1])),
  SHA  = unique(c(final_df$Taxonomie[final_df$SHA_SAR == 1],
                  final_df$Taxonomie[final_df$SHA_SCH == 1])),
  GLA  = unique(c(final_df$Taxonomie[final_df$GLA_FUG == 1],
                  final_df$Taxonomie[final_df$GLA_SME == 1],
                  final_df$Taxonomie[final_df$GLA_LIL == 1],
                  final_df$Taxonomie[final_df$GLA_KON == 1],
                  final_df$Taxonomie[final_df$GLA_DAH == 1],
                  final_df$Taxonomie[final_df$GLA_YME == 1])),
  HAR  = unique(c(final_df$Taxonomie[final_df$HAR_NYA == 1],
                  final_df$Taxonomie[final_df$HAR_BAR == 1],
                  final_df$Taxonomie[final_df$HAR_LON == 1])),
  WAL  = unique(c(final_df$Taxonomie[final_df$WAL_POO == 1],
                  final_df$Taxonomie[final_df$WAL_GUL == 1],
                  final_df$Taxonomie[final_df$WAL_SME == 1]))
)

# Function to obtain taxa that are unique to a given category
# (i.e. absent from all other categories)
get_unique_taxa <- function(category_name, taxon_lists) {
  current <- taxon_lists[[category_name]]
  others  <- unlist(taxon_lists[names(taxon_lists) != category_name])
  setdiff(current, others)
}

library(purrr)
# Compute unique taxa for each category
unique_taxa <- map(names(taxon_lists_grouped),
                   ~ get_unique_taxa(.x, taxon_lists_grouped))
names(unique_taxa) <- names(taxon_lists_grouped)

# Output: taxa exclusive to each environmental category
unique_taxa


# Intersections for all possible category combinations
# This identifies taxa shared by two or more categories
category_names <- names(taxon_lists_grouped)
intersections <- list()

for (n in 2:length(category_names)) {
  combs <- combn(category_names, n, simplify = FALSE)
  for (comb in combs) {
    shared <- Reduce(intersect, taxon_lists_grouped[comb])
    if (length(shared) > 0) {
      intersections[[paste(comb, collapse = "_AND_")]] <- shared
    }
  }
}

# Outputs
unique_taxa    # Taxa unique to each category
intersections  # Taxa shared among category combinations


#### BARPLOT ####
# Visualization of taxonomic richness per site and per vertebrate class

library(dplyr)
library(tidyr)
library(ggplot2)

# Columns corresponding to site presence/absence variables
site_cols <- 19:34

# Convert data to long format for counting
long_df <- final_df %>%
  pivot_longer(cols = all_of(site_cols), names_to = "Site", values_to = "present")

# Vertebrate classes included in the analysis
classes_of_interest <- c("Aves", "Elasmobranchii", "Mammalia", "Teleostei")

# Count the number of taxa per site and per class
count_df <- long_df %>%
  filter(present == 1, Class %in% classes_of_interest) %>%
  group_by(Site, Class) %>%
  summarise(nb_taxons = n(), .groups = "drop")

# Enforce site order according to column order in the data frame
count_df$Site <- factor(count_df$Site, levels = names(final_df)[site_cols])

# Ensure all Site × Class combinations are present (fill missing with 0)
count_df_full <- count_df %>%
  tidyr::complete(Site, Class, fill = list(nb_taxons = 0)) %>%
  mutate(Class = factor(Class, levels = c("Teleostei", "Mammalia", "Aves", "Elasmobranchii")))

# Barplot: number of taxa per site, split by class
ggplot(count_df_full, aes(x = Site, y = nb_taxons, fill = Class)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black") +
  scale_fill_brewer(palette = "Set2", name = "Class") +
  labs(
    x = "Sampling Site",
    y = "Number of Taxa"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text = element_text(color = "black"),
    legend.position = "right",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
  )


# Focus on Mammalia detections only
mammalia_taxa <- long_df %>%
  filter(Class == "Mammalia", present == 1) %>%
  mutate(
    # Prepare species names for ordered plotting
    Taxonomie = factor(Taxonomie),
    Taxonomie = forcats::fct_rev(Taxonomie)
  )

# Visual observations (red points in downstream plots)
obs_visu <- tibble::tribble(
  ~Taxonomie,                     ~Site,       ~type,
  "Odobenus rosmarus",            "WAL_POO",  "visual",
  "Odobenus rosmarus",            "SHA_SAR",  "visual",
  "Odobenus rosmarus",            "WAL_GUL",  "visual",
  "Phoca vitulina",               "SHA_SAR",  "visual",
  "Phoca vitulina",               "SCHA_SCH", "visual",
  "Phoca vitulina",               "GLA_KON",  "visual",
  "Balaenoptera acutorostrata",   "GLA_FUG",  "visual"
)

# Acoustic observations (green points in downstream plots)
obs_acou <- tibble::tribble(
  ~Taxonomie,                    ~Site,      ~type,
  "Erignathus barbatus",         "GLA_FUG", "acoustic"
)

# Harmonize factor levels for consistent plotting
obs_visu$Taxonomie <- factor(obs_visu$Taxonomie, levels = levels(mammalia_taxa$Taxonomie))
obs_acou$Taxonomie <- factor(obs_acou$Taxonomie, levels = levels(mammalia_taxa$Taxonomie))




#### MAMMALIA DETECTION HEATMAP (eDNA + VISUAL + ACOUSTIC) ####
# This figure combines eDNA detections (tiles) with independent visual and
# acoustic observations (points) for marine mammals across sampling sites.

# Heatmap of Mammalia detections using eDNA, with overlays for other methods
ggplot(mammalia_taxa, aes(x = Site, y = Taxonomie)) +
  
  ## --- eDNA detections (tiles) ---
  # Presence/absence is shown as colored tiles
  geom_tile(
    aes(fill = as.factor(present)),
    color = "grey70",
    linewidth = 0.4
  ) +
  
  scale_fill_manual(
    values = c("1" = "#4C72B0"),
    name  = "Detection type",
    labels = c("1" = "eDNA detected")
  ) +
  
  ## --- Visual observations ---
  # Independent visual sightings (red circles)
  geom_point(
    data = obs_visu,
    aes(x = Site, y = Taxonomie, shape = "Visual observation"),
    color = "red",
    size = 3.5,
    alpha = 0.9
  ) +
  
  ## --- Acoustic observations ---
  # Independent acoustic detections (green triangles)
  geom_point(
    data = obs_acou,
    aes(x = Site, y = Taxonomie, shape = "Acoustic detection"),
    color = "green3",
    size = 3.5,
    alpha = 0.9
  ) +
  
  # Shape legend for non-eDNA observations
  scale_shape_manual(
    name = "Observation type",
    values = c(
      "Visual observation" = 16,
      "Acoustic detection" = 17
    )
  ) +
  
  # Italicized species names on the y-axis
  scale_y_discrete(
    labels = function(x) parse(text = paste0("italic('", x, "')"))
  ) +
  
  labs(
    x = "Sampling site",
    y = NULL
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold"),
    axis.text.y  = element_text(face = "italic"),
    axis.title.x = element_text(face = "bold", margin = margin(t = 10)),
    panel.grid   = element_blank(),
    axis.ticks   = element_blank(),
    plot.margin  = margin(10, 10, 10, 10),
    legend.position = "right"
  )

# Export the figure at high resolution
ggsave(
  filename = "mammalia_detection_svalbard_highres.png",
  width = 10,
  height = 6,
  dpi = 600
)


#### TAXONOMIC BETA-DIVERSITY (JACCARD) ####
# This section computes taxonomic beta-diversity among sites using the
# Jaccard dissimilarity index, including decomposition into turnover and
# nestedness components.

library(textshape)
library(vegan)
library(dplyr)
library(adespatial)

# 1. Extract taxonomy and site presence/absence columns
df_sites <- final_df %>%
  select(Taxonomie, all_of(19:ncol(final_df)))

# 2. Ensure unique taxon names (required for row names)
df_sites <- df_sites %>%
  mutate(Taxonomie = make.unique(as.character(Taxonomie), sep = "_"))

# 3. Convert taxonomy column to row names
df_matrix <- df_sites %>%
  column_to_rownames("Taxonomie")

# 4. Convert to binary presence/absence matrix
jaccard <- (df_matrix > 0) * 1

# 5. Transpose so that sites are rows
jaccard_sites <- t(jaccard)

# 6. Compute Jaccard beta-diversity using adespatial
JAC <- beta.div.comp(jaccard_sites, quant = FALSE, save.abc = FALSE)
JAC

# Extract Jaccard dissimilarity matrix
JACD <- JAC$D
JACD

# Convert dist object to matrix
JACD_mat <- as.matrix(JACD)

# Optional: convert to data.frame for export
JACD_df <- as.data.frame(JACD_mat)

# Export Jaccard dissimilarity matrix to Excel
library(openxlsx)
write.xlsx(JACD_df, "JACD_matrix.xlsx", rowNames = TRUE)

# Summary statistics for pairwise dissimilarities
site_pairs <- which(upper.tri(JACD_mat), arr.ind = TRUE)
distances  <- JACD_mat[upper.tri(JACD_mat)]

min_dist  <- min(distances)
max_dist  <- max(distances)
mean_dist <- mean(distances)

# Identify site pairs with minimum and maximum dissimilarity
min_pair <- site_pairs[which.min(distances), ]
max_pair <- site_pairs[which.max(distances), ]

cat("Minimum distance :", min_dist,
    "between", rownames(JACD_mat)[min_pair[1]],
    "and", colnames(JACD_mat)[min_pair[2]], "\n")
cat("Maximum distance :", max_dist,
    "between", rownames(JACD_mat)[max_pair[1]],
    "and", colnames(JACD_mat)[max_pair[2]], "\n")
cat("Mean distance :", mean_dist, "\n")

# Decomposition into similarity, turnover, and nestedness
JAC_sim  <- 1 - JACD
JACturn <- JAC$repl
JACnest <- JAC$rich


#### PHYLOGENETIC BETA-DIVERSITY ####
# This section computes phylogenetic beta-diversity using UniFrac and
# phylogenetic Jaccard indices based on a synthesized vertebrate phylogeny.

# Build a taxonomic table for phylogeny reconstruction
taxa_list <- tibble::tibble(
  species = final_df$Taxonomie,
  genus   = final_df$Genus,
  family  = final_df$Family,
  order   = final_df$Order
)

# Prepare species list for U.PhyloMaker
sp.list <- taxa_list %>%
  dplyr::select(species, genus, family, order) %>%
  dplyr::mutate(
    species.relative = NA,
    genus.relative   = NA
  )

# Load megatrees for major vertebrate groups
library(ape)
megatree_fish  <- read.tree("phylogeny/fish_megatree.tre")
megatree_mm    <- read.tree("phylogeny/mammal_megatree.tre")
megatree_oi    <- read.tree("phylogeny/bird_megatree.tre")
megatree_shark <- read.tree("phylogeny/shark_megatree_consensus.tre")

# Combine all megatrees into a single phylogeny
combined_tree <- bind.tree(megatree_fish, megatree_mm)
combined_tree <- bind.tree(combined_tree, megatree_oi)
combined_tree <- bind.tree(combined_tree, megatree_shark)

# Load genus reference lists
gen.list_fish  <- read.csv("phylogeny/fish_genus_list.csv")
gen.list_mm    <- read.csv("phylogeny/mammal_genus_list.csv")
gen.list_oi    <- read.csv("phylogeny/bird_genus_list.csv")

# Merge genus lists
gen.list <- rbind(gen.list_fish, gen.list_mm, gen.list_oi)

# Generate a phylogeny matching the detected species
library(U.PhyloMaker)
result <- phylo.maker(sp.list, combined_tree, gen.list,
                      nodes.type = 1, scenario = 3)

phylo <- result$phylo

# Correct negative branch lengths if present
phylo$edge.length[phylo$edge.length < 0] <- 1e-6

# Root the tree using midpoint rooting
library(phangorn)
phylo_rooted <- midpoint(phylo)

# Visualize the phylogeny
library(ggtree)
ggtree(phylo_rooted, layout = "rectangular") + geom_tiplab(size = 2)
ggtree(phylo_rooted, layout = "circular")    + geom_tiplab(size = 2)

# Prepare site-by-species matrix for phylogenetic beta-diversity
jaccard_par_annee <- t(jaccard)
colnames(jaccard_par_annee) <- gsub(" ", "_", colnames(jaccard_par_annee))

# Remove selected sites prior to analysis
sites_to_remove <- c("GLA_KON", "DEEP_SFN", "GLA_YME")
jaccard_par_annee <- jaccard_par_annee[!rownames(jaccard_par_annee) %in% sites_to_remove, ]

# Compute UniFrac distances
library(picante)
unifrac <- unifrac(jaccard_par_annee, phylo_rooted)
unifrac_sim <- 1 - unifrac

# Compute phylogenetic Jaccard beta-diversity
library(betapart)
jaccard_par_annee <- jaccard_par_annee[, -1]
phylobeta <- phylo.beta.pair(jaccard_par_annee, phylo_rooted,
                             index.family = "jaccard")

# Extract and export phylogenetic Jaccard distance matrix
jac_dist <- phylobeta$phylo.beta.jac
jac_mat  <- as.matrix(jac_dist)
jac_df   <- as.data.frame(jac_mat)

library(openxlsx)
write.xlsx(jac_df, "phylo_beta_jac.xlsx", rowNames = TRUE)


# Extract and convert to matrix
jac_mat <- as.matrix(phylobeta$phylo.beta.jac)

# Select only the upper triangle to avoid duplicates
site_pairs <- which(upper.tri(jac_mat), arr.ind = TRUE)
distances <- jac_mat[upper.tri(jac_mat)]

# Calculate min, max, and mean distances
min_dist <- min(distances)
max_dist <- max(distances)
mean_dist <- mean(distances)

# Identify the corresponding site pairs
min_pair <- site_pairs[which.min(distances), ]
max_pair <- site_pairs[which.max(distances), ]

# Display results
cat("Minimum distance:", min_dist, "between", rownames(jac_mat)[min_pair[1]], "and", colnames(jac_mat)[min_pair[2]], "\n")
cat("Maximum distance:", max_dist, "between", rownames(jac_mat)[max_pair[1]], "and", colnames(jac_mat)[max_pair[2]], "\n")
cat("Mean distance:", mean_dist, "\n")

# Extract turnover and nestedness matrices
phylobetaturn <- phylobeta$phylo.beta.jtu
phylobetanest<- phylobeta$phylo.beta.jne

# Compute global sums for turnover, nestedness, and total dissimilarity
total_turnover <- sum(phylobetaturn)  # Total turnover contribution
total_nestedness <- sum(phylobetanest)  # Total nestedness contribution
total_dissimilarity <- total_turnover + total_nestedness  # Total phylogenetic dissimilarity

# Calculate percentages
percent_turnover <- (total_turnover / total_dissimilarity) * 100
percent_nestedness <- (total_nestedness / total_dissimilarity) * 100

# Display results
cat("Global percentage of turnover (phylo.beta.jtu):", round(percent_turnover, 2), "%\n")
cat("Global percentage of nestedness (phylo.beta.jne):", round(percent_nestedness, 2), "%\n")


#### PCA beta-diversity ####
#### PCA β-diversity (publication-ready) ####

library(factoextra)
library(ggplot2)

# ---- 1. Function to calculate mean pairwise β-diversity ----
mean_pairwise_beta <- function(mat) {
  mat <- as.matrix(mat)
  station_means <- rowMeans(mat, na.rm = TRUE)
  return(station_means)
}

# ---- 2. Calculate mean values for each β-diversity component ----
Bjne_Phyl <- mean_pairwise_beta(phylobetanest)
Bjtu_Phyl <- mean_pairwise_beta(phylobetaturn)
Bjtu_Tax  <- mean_pairwise_beta(JACturn)
Bjne_Tax  <- mean_pairwise_beta(JACnest)

# ---- 3. Combine into a single matrix ----
combined_matrix <- cbind(
  Bjne_Phyl = Bjne_Phyl,
  Bjtu_Phyl = Bjtu_Phyl,
  Bjtu_Tax  = Bjtu_Tax,
  Bjne_Tax  = Bjne_Tax
)

# ---- 4. Station names ----
statnam <- c("HAR_BAR", "GLA_DAH", "GLA_FUG", "WAL_GUL", "DEEP_KON", 
             "GLA_LIL", "HAR_LON", "HAR_NYA", "WAL_POO", "SHA_SAR", 
             "SHA_SCH", "WAL_SME", "GLA_SME")

# ---- 5. Filter out sites if needed ----
exclude_sites <- c("GLA_KON", "DEEP_SFN", "GLAY_YME") # adjust if necessary
statnam_filtered <- statnam[!statnam %in% exclude_sites]
combined_matrix_filtered <- combined_matrix[statnam_filtered, ]

# ---- 6. Group sites by prefix ----
groupes <- sub("_.*", "", statnam_filtered)

# ---- 7. PCA computation ----
res.pca <- prcomp(combined_matrix_filtered, scale. = TRUE)


# ---- 8. Publication-ready biplot ----
RES.PCA <- fviz_pca_biplot(
  res.pca,
  habillage = groupes,             # color points by group
  geom.ind = c("point", "text"), 
  pointshape = 21,
  pointsize = 4,
  fill.ind = groupes,              # fill color for points
  col.ind = groupes,               # text color (labels)
  mean.point = FALSE,              # do not show group centroids
  label = "all",                   # display all station labels
  repel = TRUE,                    # avoid overlapping labels
  col.var = "black",               # variable arrows in black
  arrowsize = 0.8,
  labelsize = 4,
  legend.title = "Habitat",
  title = NULL,
  addEllipses = TRUE,
  ellipse.type = "convex"
) +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Set2") +  # same palette for points and text
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14, color = "black"),
    panel.grid = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )
RES.PCA
# ---- 9. Save the plot ----
ggsave(
  filename = "PCA_beta_diversity_publication_ready.png",
  plot = RES.PCA,
  width = 14,
  height = 8,
  dpi = 600
)


#### Phylogenetic Diversity (PD) calculations ####
#### PD metrics by categories ####

# 1. Extract row names
site_names <- rownames(jaccard_par_annee)

# 2. Extract categories (GLA, WAL, etc.)
categories <- sub("_.*", "", site_names)

# 3. Convert to data.frame and add category column
df <- as.data.frame(jaccard_par_annee)
df$category <- categories

# 4. Group by category and take max or mean (depending on need)
library(dplyr)

# Example here with max (presence/absence in the category)
df_grouped <- df %>%
  group_by(category) %>%
  summarise(across(.cols = everything(), .fns = max)) %>%
  as.data.frame()
df_grouped

# 5. Binarize: >0 -> 1, else 0
df_binary <- as.data.frame((df_grouped[,-1] > 0) * 1)

# 6. Restore row names (categories)
rownames(df_binary) <- df_grouped$category

# Harmonize column names
colnames(df_binary) <- gsub(" ", "_", colnames(df_binary))

# Intersect species present in matrix and phylogeny
species_common <- intersect(colnames(df_binary), phylo_rooted$tip.label)

jaccard_filt <- df_binary[, species_common]
phydist_filt <- cophenetic(phylo_rooted)[species_common, species_common]

phylo <- result$phylo
phylo$edge.length[phylo$edge.length < 0] <- 1e-6
phylo_rooted <- midpoint(phylo)

#### PD / SES calculations ####
pd.result <- pd(jaccard_filt, phylo_rooted, include.root = TRUE)

ses.pd.result <- picante::ses.pd(jaccard_filt, phylo, null.model = "taxa.labels", runs = 1000, include.root=FALSE)
ses.pd.result
ses.mpd.result <- ses.mpd(jaccard_filt, phydist_filt, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 1000)
ses.mntd.result <- ses.mntd(jaccard_filt, phydist_filt, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 1000)

#### Merge results ####
# Observed values
recap_obs <- pd.result %>%
  as.data.frame() %>%
  mutate(Site = rownames(.)) %>%
  left_join(
    ses.mpd.result %>% as.data.frame() %>% mutate(Site = rownames(.)) %>% select(Site, mpd.obs),
    by = "Site"
  ) %>%
  left_join(
    ses.mntd.result %>% as.data.frame() %>% mutate(Site = rownames(.)) %>% select(Site, mntd.obs),
    by = "Site"
  ) %>%
  select(Site, SR, PD, mpd.obs, mntd.obs)

# Long format for plotting Panel A
recap_obs_long <- recap_obs %>%
  pivot_longer(cols = c(SR, PD, mpd.obs),
               names_to = "Metric", values_to = "Observed") %>%
  mutate(Metric = recode(Metric, "SR"="SR", "PD"="PD", "mpd.obs"="MPD"),
         Metric = factor(Metric, levels = c("SR","PD","MPD")))

# SES results (z-score + p-value)
ses_long <- bind_rows(
  ses.pd.result %>% as.data.frame() %>% mutate(Site=rownames(.)) %>% select(Site, z=pd.obs.z, p=pd.obs.p) %>% mutate(Metric="PD"),
  ses.mpd.result %>% as.data.frame() %>% mutate(Site=rownames(.)) %>% select(Site, z=mpd.obs.z, p=mpd.obs.p) %>% mutate(Metric="MPD"),
) %>%
  mutate(Significant = ifelse(p < 0.05 | p > 0.95, "Yes", "No"))
ses_long


#### Plots ####
# Panel A: Observed metrics
plot_obs <- ggplot(recap_obs_long, aes(x=Site, y=Observed, fill=Metric)) +
  geom_col() +
  facet_wrap(~ Metric, scales="free_x") +
  coord_flip() +
  theme_classic(base_size=14) +
  labs(title="Observed phylogenetic metrics", x="Site", y="Observed value") +
  scale_fill_brewer(palette="Set2") +
  theme(legend.position="none")

# Panel B: SES z-scores
plot_ses <- ggplot(ses_long, aes(x=reorder(Site, z), y=z, color=Metric)) +
  geom_hline(yintercept=c(-1.96,1.96), linetype="dashed", color="grey40") +
  geom_point(aes(shape=Significant), size=3, position=position_dodge(width=0.6)) +
  coord_flip() +
  theme_classic(base_size=14) +
  labs(title="Standardized Effect Sizes (SES)", x="Site", y="SES (z-score)") +
  scale_color_brewer(palette="Set2") +
  scale_shape_manual(values=c("Yes"=17,"No"=16)) +
  theme(legend.position="top")


# Combine panels
library(patchwork)
final_plot <- plot_obs + plot_ses + plot_layout(ncol=2, widths=c(1.3,1.6)) + plot_annotation(tag_levels="A")
final_plot


library(dplyr)
library(fmsb)

# Prepare data: SES z-scores + Species Richness (SR)
ses_for_plot <- ses_long %>%
  select(Site, Metric, z) %>%
  pivot_wider(names_from = Metric, values_from = z) %>%
  left_join(recap_obs %>% select(Site, SR), by="Site") %>%
  rename(PD = PD, MPD = MPD)  # PD here can be SES or observed, adjust if needed

# Replace NA values with 0
radar_data <- ses_for_plot %>% 
  mutate(across(c(SR, PD, MPD), ~replace_na(., 0)))

# Normalize each metric between 0 and 1
radar_data_scaled <- radar_data %>%
  mutate(across(c(SR, PD, MPD), ~(. - min(., na.rm=TRUE)) / (max(., na.rm=TRUE) - min(., na.rm=TRUE))))

# Define max and min for fmsb radar chart
max_min <- data.frame(
  SR = c(1,0),
  PD = c(1,0),
  MPD = c(1,0)
)

# Combine max/min and scaled data
radar_plot_data <- rbind(max_min, radar_data_scaled %>% select(SR, PD, MPD))
rownames(radar_plot_data) <- c("max","min", radar_data_scaled$Site)

# Save individual radar charts as PDF
pdf("radar_sites_scaled.pdf", width=10, height=7)
par(mar = c(1,1,2,1))
for(i in 3:nrow(radar_plot_data)){
  radarchart(radar_plot_data[c(1,2,i),],
             axistype=1,
             pcol="blue", pfcol=rgb(0,0,1,0.3), plwd=2,
             cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,1,0.2), cglwd=0.8,
             vlcex=0.8)
  mtext(rownames(radar_plot_data)[i], side=3, line=1, cex=1.2, font=2)  # Add site name above radar
}
dev.off()


library(fmsb)

# Number of sites (excluding max/min rows)
n_sites <- nrow(radar_plot_data) - 2
sites <- rownames(radar_plot_data)[-c(1,2)]

# Define grid layout: 3 columns, calculate number of rows
n_col <- 3
n_row <- ceiling(n_sites / n_col)

# Create PDF for grid of radar charts
pdf("radar_sites_grid.pdf", width = 14, height = 8)
par(mfrow = c(n_row, captucapcan_col), mar = c(1,1,2,1))  # Smaller margins for mini-radars

# Loop to plot radar for each site
for(i in 1:n_sites){
  radarchart(radar_plot_data[c(1,2,i+2),],  # +2 to skip max/min rows
             axistype=2,
             pcol="blue", pfcol=rgb(0,0,1,0.3), plwd=2,
             cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,1,0.2), cglwd=0.8,
             vlcex=1.5)
  mtext(sites[i], side=3, line=0.5, cex=0.9, font=2)  # Add site name above radar
}

dev.off()


library(fmsb)

# 1️⃣ Select 5 key sites
sites_to_plot <- c("GLA", "HAR", "WAL", "DEEP", "SHA")
radar_plot_data_subset <- radar_plot_data[c("max","min", sites_to_plot), ]

# 2️⃣ Variable names for radar axes
variable_names <- c("TR", "SES.PD", "SES.MPD")

# ✅ 2b. Construct general radar chart (all values = 1)
radar_general <- rbind(
  max = rep(1, length(variable_names)),
  min = rep(0, length(variable_names)),
  general = rep(1, length(variable_names))
)
colnames(radar_general) <- variable_names
radar_general <- as.data.frame(radar_general)

# 3️⃣ Save PNG
png("radar_plots.png", width = 2000, height = 1500, res = 200)

# 4️⃣ Grid layout: 2 rows x 3 columns
par(mfrow = c(2,3), mar = c(2,1,4,2))

# 5️⃣ Plot radar charts for individual sites
for(i in 3:nrow(radar_plot_data_subset)){
  radarchart(
    radar_plot_data_subset[c(1,2,i), ],
    axistype = 1,
    pcol = "blue",
    pfcol = rgb(0,0,1,0.3),
    plwd = 2,
    cglcol = "grey",
    cglty = 1,
    axislabcol = "grey",
    caxislabels = seq(0,1,1),
    cglwd = 0.8,
    vlcex = 1.2,
    vlabels = variable_names
  )
  
  mtext(
    rownames(radar_plot_data_subset)[i], 
    side = 3, 
    line = 1,
    cex = 1.5,
    font = 2
  )
}

# ✅ 6️⃣ Plot the general radar chart (6th panel)
radarchart(
  radar_general,
  axistype = 1,
  pcol = "darkred",
  pfcol = rgb(1,0,0,0.3),
  plwd = 3,
  cglcol = "grey",
  cglty = 1,
  axislabcol = "grey",
  caxislabels = seq(0,1,0.2),
  cglwd = 0.8,
  vlcex = 1.2
)

mtext("General (all values = 1)", side = 3, line = 1, cex = 1.5, font = 2)

# 7️⃣ Reset layout and close PNG device
par(mfrow = c(1,1))
dev.off()



#### PHYLOGENETIC TREE – PUBLICATION VERSION ####

library(ggtree)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(stringr)

# --- 1. Prepare data ---
site_cols <- final_df %>%
  select(where(is.integer)) %>%
  colnames()

data_long <- final_df %>%
  select(Taxonomie, Class, all_of(site_cols)) %>%
  pivot_longer(
    cols = all_of(site_cols),
    names_to = "site",
    values_to = "presence"
  ) %>%
  filter(presence > 0) %>%
  mutate(
    site_grouped = str_extract(site, "^[^_]+"),  # Extract site prefix
    Taxonomie = str_replace_all(Taxonomie, " ", "_")  # Replace spaces with underscore
  ) %>%
  select(Taxonomie, Class, site_grouped)

sites <- unique(data_long$site_grouped)

# --- 2. Define colors per site ---
station_colors <- c(
  "HAR"  = "blue",
  "GLA"  = "lightblue",
  "WAL"  = "orange",
  "DEEP" = "purple",
  "SHA"  = "green"
)

# --- 3. Calculate positions of orders for clade annotation ---
if(!"order" %in% colnames(data_long)) stop("Add 'order' column to data_long")

ord_positions <- lapply(unique(data_long$order), function(ord) {
  tips_ord <- data_long %>% filter(order == ord) %>% pull(Taxonomie)
  tips_ord <- tips_ord[tips_ord %in% phylo_rooted$tip.label]
  if(length(tips_ord) > 1) {
    mrca <- getMRCA(phylo_rooted, tips_ord)  # Most recent common ancestor
    if(!is.null(mrca)) return(data.frame(order = ord, node = mrca))
  }
  NULL
}) %>% bind_rows()

# --- 4. Function to plot a phylogenetic tree per site ---
plot_phylo_site <- function(site_name, color_map) {
  
  tip_df <- data.frame(label = phylo_rooted$tip.label) %>%
    mutate(present_in_site = ifelse(
      label %in% data_long$Taxonomie[data_long$site_grouped == site_name],
      site_name, NA
    ))
  
  p <- suppressWarnings(
    ggtree(phylo_rooted, layout = "rectangular", ladderize = FALSE) %<+% tip_df +
      geom_tippoint(
        aes(color = present_in_site),
        size = 1.5,
        na.rm = TRUE
      ) +
      scale_color_manual(values = color_map, na.value = NA, name = "Site") +
      theme_void() +
      ggtitle(site_name) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        legend.position = "none",
        plot.margin = margin(10, 10, 10, 10)
      ) +
      coord_cartesian(clip = "off")
  )
  
  # Add clade labels
  offset_base <- 3.5
  offset_step <- ifelse(nrow(ord_positions) > 1, 3 / nrow(ord_positions), 0.3)
  
  for(i in seq_len(nrow(ord_positions))) {
    p <- p + geom_cladelab(
      node = ord_positions$node[i],
      label = ord_positions$order[i],
      align = TRUE,
      offset = offset_base + (i - 1) * offset_step,
      barsize = 0.8,
      fontsize = 3.5,
      hjust = 0
    )
  }
  
  # --- Optional: rotate node 36 ---
  p_rot <- tryCatch(rotate(p, node = 36), error = function(e) p)
  
  return(p_rot)
}


# --- 5. Generate tree for each site ---
plots_sites <- lapply(sites, plot_phylo_site, color_map = station_colors)
plots_sites

# --- 6. Create full tree with same orientation ---
tip_df_all <- data.frame(
  label = phylo_rooted$tip.label,
  present_in_site = "All Taxa"
)

full_tree <- suppressWarnings(
  ggtree(phylo_rooted, layout = "rectangular", ladderize = FALSE) %<+% tip_df_all +
    geom_tiplab(size = 2.3) +       # show all tip labels
    theme_void() +
    ggtitle("All Taxa") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      plot.margin = margin(10, 10, 10, 10)
    ) +
    coord_cartesian(clip = "off")
)

full_tree <- suppressWarnings(
  ggtree(phylo_rooted, layout = "rectangular", ladderize = FALSE) %<+% tip_df_all +
    geom_tiplab(size = 2.3) +
    expand_limits(x = max(phylo_rooted$edge.length) * 1.5) +   # extra space
    theme_void() +
    ggtitle("All Taxa") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      plot.margin = margin(10, 40, 10, 10)  # increase right margin
    ) +
    coord_cartesian(clip = "off")
)

full_tree_swapped <- ggtree::rotate(full_tree, node = 36)
full_tree_swapped

# --- 7. Combine all site plots ---
plots_sites_combined <- wrap_plots(plots_sites, ncol=5)
plots_sites_combined

# Combine site-specific trees with full tree
combined_plot <- (plots_sites_combined/full_tree_swapped)
combined_plot

# --- 8. Save and display ---
ggsave(
  filename = "Figure_phylo_tree_by_site_PolarBiology.png",
  plot = combined_plot,
  width = 10, height = 10, dpi = 300
)

print(combined_plot)

