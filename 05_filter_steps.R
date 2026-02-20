
load_packages <- function(pkg_list) {
  to_install <- pkg_list[!pkg_list %in% installed.packages()[, "Package"]]
  if (length(to_install)) install.packages(to_install, dependencies = TRUE)
  invisible(lapply(pkg_list, library, character.only = TRUE))
}

# List of packages commonly used for QIIME2 + Phyloseq analysis
packages <- c(
  "qiime2R",      # Import QIIME2 artifacts (.qza/.qzv)
  "phyloseq",     # Microbiome data handling
  "tidyverse",    # Data manipulation (dplyr, ggplot2, etc.)
  "vegan",        # Diversity metrics and ordination
  "microbiome",   # Extra microbiome analysis tools
  "RColorBrewer", # Color palettes
  "plotly",        # Interactive plots (optional)
  "patchwork"
)

# Load (and install if missing)
load_packages(packages)

PS_16S<-qza_to_phyloseq(
  features="C:\\Users\\mhewapat\\Documents\\avery_Farms\\16S_full_run\\results\\artifacts\\table_16S.qza",
  taxonomy="C:\\Users\\mhewapat\\Documents\\avery_Farms\\16S_full_run\\results\\artifacts\\taxonomy_16S_silva138_99.qza",
  metadata="C:\\Users\\mhewapat\\Documents\\avery_Farms\\16S_full_run\\results\\inputs\\metadata_clean.tsv"
)

PS_16S
nsamples(PS_16S)
ntaxa(PS_16S)

sample_reads <- sample_sums(PS_16S)

summary(sample_reads)
quantile(sample_reads, probs = c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1))

asv_total_reads <- taxa_sums(PS_16S)
total_reads <- sum(asv_total_reads)

# Convert to percent of total reads
asv_percent_total <- (asv_total_reads / total_reads) * 100

summary(asv_percent_total)

sum(asv_percent_total < 0.01)
sum(asv_percent_total < 0.005)
sum(asv_percent_total < 0.1)

length(asv_percent_total)

hist(log10(asv_total_reads + 1),
     breaks = 50,
     main = "ASV abundance distribution (log10 scale)",
     xlab = "log10(total ASV reads + 1)")

sum(asv_percent_total < 0.01)

table(asv_prevalence)[1:5]

asv_prevalence <- apply(otu_table(PS_16S), 1, function(x) sum(x > 0))

summary(asv_prevalence)
table(asv_prevalence)[1:10]   # how many ASVs appear in only 1, 2, 3 samples

n_samp <- nsamples(PS_16S)
asv_prev_percent <- (asv_prevalence / n_samp) * 100

summary(asv_prev_percent)

library(phyloseq)

## ---------------------------
## 1. Starting stats
## ---------------------------
start_asvs  <- ntaxa(PS_16S)
start_reads <- sum(sample_sums(PS_16S))

## ---------------------------
## 2. Global abundance filter (0.01%)
## ---------------------------
asv_total_reads <- taxa_sums(PS_16S)
total_reads     <- sum(asv_total_reads)

total_reads

asv_percent_total <- (asv_total_reads / total_reads) * 100

asv_percent_total

PS_abund_filt <- prune_taxa(asv_percent_total >= 0.001, PS_16S)

## ---------------------------
## 3. Summary
## ---------------------------
summary_df <- data.frame(
  Metric = c("ASVs", "Samples", "Total Reads"),
  Before = c(start_asvs, nsamples(PS_16S), start_reads),
  After  = c(ntaxa(PS_abund_filt), nsamples(PS_abund_filt), sum(sample_sums(PS_abund_filt)))
)

summary_df

(sum(sample_sums(PS_abund_filt)) / start_reads) * 100

###############################################################################
# Alpha diversity (Observed ASVs) comparison: Leaf vs Sponge
#
# Purpose:
#   - Compare species richness (Observed ASVs) between Leaf and Sponge samples
#   - Uses rarefied data to control for sequencing depth
#   - Performs a Wilcoxon rank-sum test (non-parametric)
#
# Input:
#   - PS_16S_Sponge_no_control (phyloseq object, no controls)
#
# Output:
#   - Boxplot with individual samples overlaid
#   - Wilcoxon p-value comparing Leaf vs Sponge
###############################################################################

library(phyloseq)
library(tidyverse)
library(ggpubr)

# 1) Rarefy to minimum sequencing depth
set.seed(123)

PS_abund_filt

PS_16S

PS_rare <- rarefy_even_depth(
  PS_abund_filt,
  rngseed = 123,
  sample.size = 29600,
  replace = FALSE,
  verbose = FALSE
)

sample_sums(PS_rare)

# 2) Calculate Observed ASVs
alpha_df <- estimate_richness(
  PS_rare,
  measures = "Observed"
)

# 3) Add metadata
alpha_df <- estimate_richness(PS_rare, measures = "Observed") %>%
  tibble::rownames_to_column("Sample")

meta_df <- data.frame(sample_data(PS_rare)) %>%
  tibble::rownames_to_column("Sample")

alpha_df <- left_join(alpha_df, meta_df, by = "Sample")


# 4) Plot + statistical test
ggplot(alpha_df, aes(x = Sample_Type, y = Observed, fill = Sample_Type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = list(c("Leaf", "Sponge"))
  ) +
  scale_fill_manual(values = c("Leaf" = "#028A0F", "Sponge" = "#02A3D3"))
theme_bw(base_size = 13) +
  labs(
    title = "Observed ASVs by sample type",
    x = "Sample Type",
    y = "Observed ASVs"
  ) +
  theme(legend.position = "none")


###############################################################################
# Beta diversity (Jaccard presence/absence) — Timepoint comparison
#
# Purpose:
#   - Assess differences in community composition across timepoints
#   - Uses Jaccard distance (presence/absence) on rarefied data
#   - Visualizes ordination with 95% confidence ellipses
#   - Tests differences using PERMANOVA
#
# Input:
#   - ps_sponge_rarefied_no_control (already rarefied)
###############################################################################

library(phyloseq)
library(tidyverse)
library(vegan)

# 1) Convert to presence/absence
ps_pa <- transform_sample_counts(
  PS_rare,
  function(x) as.numeric(x > 0)
)

# 2) Compute Jaccard distance
jaccard_dist <- phyloseq::distance(
  PS_rare,
  method = "jaccard"
)

# 3) PCoA ordination
ord_jaccard <- ordinate(
  PS_rare,
  method = "PCoA",
  distance = jaccard_dist
)

# 4) Extract ordination scores + metadata
ord_df <- plot_ordination(
  PS_rare,
  ord_jaccard,
  justDF = TRUE
)

colnames(ord_df)

library(ggrepel)


# 5) PCoA plot with confidence ellipses
p <- ggplot(
  ord_df,
  aes(x = Axis.1, y = Axis.2, color = Timepoint)
) +
  geom_point(size = 3, alpha = 0.9) +
  geom_text_repel(
    aes(label = SampleID),   # change if column name differs
    size = 3,
    max.overlaps = 20
  ) +
  stat_ellipse(
    aes(group = Timepoint),
    type = "t",
    level = 0.95,
    linewidth = 1
  ) +
  theme_bw(base_size = 13) +
  labs(
    title = "Jaccard beta diversity (presence/absence)",
    x = paste0("PCoA1 (", round(ord_jaccard$values$Relative_eig[1] * 100, 1), "%)"),
    y = paste0("PCoA2 (", round(ord_jaccard$values$Relative_eig[2] * 100, 1), "%)"),
    color = "Timepoint"
  )

p

meta <- sample_data(PS_abund_filt)


##############################################################################
# Beta diversity (Jaccard presence/absence) — Tank_Type comparison
###############################################################################

library(phyloseq)
library(tidyverse)
library(vegan)

# Presence / absence transform
ps_pa <- transform_sample_counts(
  PS_rare,
  function(x) as.numeric(x > 0)
)

# Jaccard distance
jaccard_dist <- phyloseq::distance(
  PS_rare,
  method = "jaccard"
)

# PCoA ordination
ord_jaccard <- ordinate(
  PS_rare,
  method = "PCoA",
  distance = jaccard_dist
)

# Extract scores + metadata
ord_df <- plot_ordination(
  PS_rare,
  ord_jaccard,
  justDF = TRUE
)

# Plot
ggplot(
  ord_df,
  aes(x = Axis.1, y = Axis.2, color = Tank_Type)
) +
  geom_point(size = 3, alpha = 0.9) +
  stat_ellipse(
    aes(group = Tank_Type),
    type = "t",
    level = 0.95,
    linewidth = 1
  ) +
  theme_bw(base_size = 13) +
  labs(
    title = "Jaccard beta diversity by tank type",
    x = paste0(
      "PCoA1 (",
      round(ord_jaccard$values$Relative_eig[1] * 100, 1),
      "%)"
    ),
    y = paste0(
      "PCoA2 (",
      round(ord_jaccard$values$Relative_eig[2] * 100, 1),
      "%)"
    ),
    color = "Tank Type"
  )



# Genera of interest
genera_of_interest <- c(
  "Azospirillum",
  "Pseudomonas",
  "Bacillus",
  "Paenibacillus",
  "Pantoea",
  "Enterobacter",
  "Aeromonas",
  "Klebsiella",
  "Ralstonia"
)

# 1) Transform to relative abundance (per sample)
ps_rel <- transform_sample_counts(
  PS_abund_filt,
  function(x) x / sum(x)
)

# 2) Melt to long format
rel_df <- psmelt(ps_rel)

# 3) Handle unclassified genus
rel_df <- rel_df %>%
  mutate(Genus = ifelse(is.na(Genus), "Unclassified", Genus))

# 4) Define plotting groups
rel_df <- rel_df %>%
  mutate(
    Genus_plot = ifelse(
      Genus %in% genera_of_interest,
      Genus,
      "Other"
    )
  )

# 5) Aggregate
rel_df_sum <- rel_df %>%
  group_by(Timepoint, Tank_Type, Genus_plot) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

# 6) Re-normalize within each Timepoint × Tank_Type
rel_df_sum <- rel_df_sum %>%
  group_by(Timepoint, Tank_Type) %>%
  mutate(Abundance = Abundance / sum(Abundance)) %>%
  ungroup()

# 7) Lock levels and colors
genus_levels <- c(genera_of_interest, "Other")

rel_df_sum$Genus_plot <- factor(rel_df_sum$Genus_plot, levels = genus_levels)

genus_colors <- c(
  setNames(
    createPalette(length(genera_of_interest), c("#000000", "#FFFFFF")),
    genera_of_interest
  ),
  "Other" = "grey80"
)

# 8) Pie plotting function
plot_timepoint_pie <- function(tp, title_text, tanks) {
  ggplot(
    filter(rel_df_sum, Timepoint == tp, Tank_Type %in% tanks),
    aes(x = "", y = Abundance, fill = Genus_plot)
  ) +
    geom_col(width = 1, color = "black", linewidth = 0.2) +
    coord_polar(theta = "y") +
    facet_wrap(~ Tank_Type) +
    scale_fill_manual(values = genus_colors) +
    theme_void(base_size = 13) +
    labs(
      title = title_text,
      fill = "Genus"
    )
}

# 9) Generate pie charts
p_pre <- plot_timepoint_pie(
  "Pre_inoculation",
  "Pre-inoculation",
  c("Control", "Inoculated")
)

p_post <- plot_timepoint_pie(
  "1_week_post_inoculation",
  "Post-inoculation",
  c("Control", "Inoculated")
)

p_harvest <- plot_timepoint_pie(
  "Harvest",
  "Harvest",
  c("Control", "Inoculated", "Existing_tank")
)

# Display plots
p_pre
p_post
p_harvest

(p_pre + p_post + p_harvest) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")



###############################################################################
# Bar charts of selected genera (plus Other) across timepoints
#
# Purpose:
#   - Compare relative abundance of genera of interest between treatments
#   - Display each genus explicitly on the x-axis (including "Other")
#   - Use non-stacked bars for clear comparison
#   - Maintain identical coloring across all figures
#
# Input:
#   - rel_df_sum (from genera-of-interest workflow)
#   - genus_colors (predefined, reused palette)
###############################################################################

library(tidyverse)

# Ensure factor order is consistent
rel_df_sum$Genus_plot <- factor(
  rel_df_sum$Genus_plot,
  levels = names(genus_colors)
)

# Plotting function
plot_timepoint_bar <- function(tp, title_text, tanks) {
  ggplot(
    filter(rel_df_sum, Timepoint == tp, Tank_Type %in% tanks),
    aes(x = Genus_plot, y = Abundance, fill = Genus_plot)
  ) +
    geom_col(
      position = position_dodge(width = 0.8),
      color = "black",
      linewidth = 0.2
    ) +
    facet_wrap(~ Tank_Type) +
    scale_fill_manual(values = genus_colors) +
    theme_bw(base_size = 13) +
    labs(
      title = title_text,
      x = "Genus",
      y = "Relative abundance"
    ) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

# Generate plots
p_pre <- plot_timepoint_bar(
  "Pre_inoculation",
  "Pre-inoculation",
  c("Control", "Inoculated")
)

p_post <- plot_timepoint_bar(
  "1_week_post_inoculation",
  "Post-inoculation",
  c("Control", "Inoculated")
)

p_harvest <- plot_timepoint_bar(
  "Harvest",
  "Harvest",
  c("Control", "Inoculated", "Existing_tank")
)

# Display
p_pre
p_post
p_harvest

(p_pre + p_post + p_harvest) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")




###############################################################################
# Non-stacked bar charts of selected genera (plus Other) across timepoints
#
# - NON-rarefied data
# - Genera of interest shown explicitly
# - All other genera collapsed into "Other"
# - Percent labels above bars
# - Fixed y-axis (0–100%) across all plots
# - Consistent coloring across figures
###############################################################################

library(phyloseq)
library(tidyverse)
library(Polychrome)
library(scales)

# ----------------------------
# Genera of interest
# ----------------------------
genera_of_interest <- c(
  "Azospirillum",
  "Pseudomonas",
  "Bacillus",
  "Paenibacillus",
  "Pantoea",
  "Enterobacter",
  "Aeromonas",
  "Klebsiella",
  "Ralstonia"
)

# ----------------------------
# Relative abundance transform
# ----------------------------
ps_rel <- transform_sample_counts(
  Ps_sponge_no_euk_no_control,
  function(x) x / sum(x)
)

# ----------------------------
# Melt + clean taxonomy
# ----------------------------
rel_df <- psmelt(ps_rel) %>%
  mutate(
    Genus = ifelse(is.na(Genus), "Unclassified", Genus),
    Genus_plot = ifelse(Genus %in% genera_of_interest, Genus, "Other")
  )

# ----------------------------
# Aggregate + renormalize
# ----------------------------
rel_df_sum <- rel_df %>%
  group_by(Timepoint, Tank_Type, Genus_plot) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  group_by(Timepoint, Tank_Type) %>%
  mutate(Abundance = Abundance / sum(Abundance)) %>%
  ungroup()

# ----------------------------
# Factor order + colors
# ----------------------------
genus_levels <- c(genera_of_interest, "Other")
rel_df_sum$Genus_plot <- factor(rel_df_sum$Genus_plot, levels = genus_levels)

genus_colors <- c(
  setNames(
    createPalette(length(genera_of_interest), c("#000000", "#FFFFFF")),
    genera_of_interest
  ),
  "Other" = "grey80"
)

# Explicitly recolor Pseudomonas (avoid grey conflict)
genus_colors["Pseudomonas"] <- "#1F78B4"

# ----------------------------
# Plotting function
# ----------------------------
plot_timepoint_bar <- function(tp, title_text, tanks) {
  ggplot(
    filter(rel_df_sum, Timepoint == tp, Tank_Type %in% tanks),
    aes(x = Genus_plot, y = Abundance, fill = Genus_plot)
  ) +
    geom_col(
      position = position_dodge(width = 0.8),
      color = "black",
      linewidth = 0.2
    ) +
    geom_text(
      aes(label = percent(Abundance, accuracy = 1)),
      vjust = -0.3,
      size = 3
    ) +
    facet_wrap(~ Tank_Type) +
    scale_fill_manual(values = genus_colors) +
    scale_y_continuous(
      limits = c(0, 1),
      labels = scales::percent_format(accuracy = 1)
    ) +
    theme_bw(base_size = 13) +
    labs(
      title = title_text,
      x = "Genus",
      y = "Relative abundance (%)"
    ) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

# ----------------------------
# Generate plots
# ----------------------------
p_pre <- plot_timepoint_bar(
  "Pre_inoculation",
  "Pre-inoculation",
  c("Control", "Inoculated")
)

p_post <- plot_timepoint_bar(
  "1_week_post_inoculation",
  "Post-inoculation",
  c("Control", "Inoculated")
)

p_harvest <- plot_timepoint_bar(
  "Harvest",
  "Harvest",
  c("Control", "Inoculated", "Existing_tank")
)

# Display
p_pre
p_post
p_harvest


(p_pre + p_post + p_harvest) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")




###############################################################################
# Beta diversity (Bray–Curtis abundance) — Timepoint comparison
#
# Purpose:
#   - Assess differences in community composition across timepoints
#   - Uses Bray–Curtis distance (abundance-weighted)
#   - Does NOT require rarefaction
#   - Visualizes ordination with 95% confidence ellipses
#   - Tests differences using PERMANOVA
#
# Input:
#   - PS_abund_filt  (or PS_16S if you decide not to abundance-filter)
###############################################################################

library(phyloseq)
library(tidyverse)
library(vegan)
library(ggrepel)

# 1) Compute Bray–Curtis distance (abundance-weighted)
bray_dist <- phyloseq::distance(
  PS_abund_filt_sponge,
  method = "bray"
)

# 2) PCoA ordination
ord_bray <- ordinate(
  PS_abund_filt_sponge,
  method = "PCoA",
  distance = bray_dist
)

# 3) Extract ordination scores + metadata
ord_df <- plot_ordination(
  PS_abund_filt_sponge,
  ord_bray,
  justDF = TRUE
)

ord_df$SampleID <- rownames(ord_df)

# 4) PCoA plot with confidence ellipses + labels
ggplot(ord_df, aes(x = Axis.1, y = Axis.2, color = Timepoint)) +
  geom_point(size = 3) +
  theme_bw() +
  ggplot(ord_df, aes(x = Axis.1, y = Axis.2, color = Timepoint)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = SampleID), size = 3) +
  theme_bw()

p_bray

PS_abund_filt_sponge

library(phyloseq)
library(ggplot2)
library(ggrepel)

# 1) Compute Bray–Curtis distance (no rarefaction needed)
bray_dist <- phyloseq::distance(
  PS_abund_filt_sponge,          # or PS_abund_filt if you choose to filter
  method = "bray"
)

# 2) PCoA ordination
ord_bray <- ordinate(
  PS_abund_filt_sponge,
  method = "PCoA",
  distance = bray_dist
)

# 3) Extract ordination scores + metadata
ord_df <- plot_ordination(
  PS_abund_filt_sponge,
  ord_bray,
  justDF = TRUE
)

# 4) Add SampleID (stored as rownames)
ord_df$SampleID <- rownames(ord_df)

ord_df

while (!is.null(dev.list())) dev.off()


# 5) Bray–Curtis PCoA plot (points + labels only)
p_bray <- ggplot(
  ord_df,
  aes(x = Axis.1, y = Axis.2, color = Timepoint)
) +
  geom_point(size = 3, alpha = 0.9) +
  geom_text_repel(
    aes(label = SampleID),
    size = 3,
    max.overlaps = 20
  ) +
  theme_bw(base_size = 13) +
  labs(
    title = "Bray–Curtis beta diversity (abundance-weighted)",
    x = paste0(
      "PCoA1 (",
      round(ord_bray$values$Relative_eig[1] * 100, 1),
      "%)"
    ),
    y = paste0(
      "PCoA2 (",
      round(ord_bray$values$Relative_eig[2] * 100, 1),
      "%)"
    ),
    color = "Timepoint"
  )

p_bray

PS_abund_filt_sponge <- subset_samples(PS_rare, Sample_Type == "Sponge")


library(phyloseq)
library(ggplot2)

bray_dist <- phyloseq::distance(PS_abund_filt_sponge, method = "bray")

ord_bray <- ordinate(
  PS_abund_filt,
  method = "PCoA",
  distance = bray_dist
)

ord_df <- plot_ordination(
  PS_16S,
  ord_bray,
  justDF = TRUE
)

ggplot(ord_df, aes(x = Axis.1, y = Axis.2, color = Timepoint)) +
  geom_point(size = 3) +
  theme_bw()

ord_df$SampleID <- rownames(ord_df)

ggplot(ord_df, aes(x = Axis.1, y = Axis.2, color = Timepoint)) +
  geom_point(size = 3) +
  geom_text(
    aes(label = SampleID),
    size = 3,
    vjust = -1
  ) +
  theme_bw()



library(phyloseq)
library(ggplot2)

# 1) Convert to presence/absence
ps_pa <- transform_sample_counts(
  PS_abund_filt_sponge,
  function(x) as.numeric(x > 0)
)

# 2) Compute Jaccard distance
jaccard_dist <- phyloseq::distance(
  PS_abund_filt_sponge,
  method = "jaccard"
)

# 3) PCoA ordination
ord_jaccard <- ordinate(
  PS_abund_filt_sponge,
  method = "PCoA",
  distance = jaccard_dist
)

# 4) Extract ordination scores + metadata
ord_df <- plot_ordination(
  PS_abund_filt_sponge,
  ord_jaccard,
  justDF = TRUE
)

# 5) Basic Jaccard PCoA plot
p_jaccard <- ggplot(
  ord_df,
  aes(x = Axis.1, y = Axis.2, color = Timepoint)
) +
  geom_point(size = 3, alpha = 0.9) +
  theme_bw(base_size = 13) +
  labs(
    title = "Jaccard beta diversity (presence/absence)",
    x = paste0(
      "PCoA1 (",
      round(ord_jaccard$values$Relative_eig[1] * 100, 1),
      "%)"
    ),
    y = paste0(
      "PCoA2 (",
      round(ord_jaccard$values$Relative_eig[2] * 100, 1),
      "%)"
    ),
    color = "Timepoint"
  )

p_jaccard


library(phyloseq)
library(ggplot2)

# 1) Compute Bray–Curtis distance
bray_dist <- phyloseq::distance(
  PS_abund_filt_sponge,
  method = "bray"
)

# 2) PCoA ordination
ord_bray <- ordinate(
  PS_abund_filt_sponge,
  method = "PCoA",
  distance = bray_dist
)

# 3) Extract ordination scores + metadata
ord_df <- plot_ordination(
  PS_abund_filt_sponge,
  ord_bray,
  justDF = TRUE
)

# 4) Bray–Curtis PCoA plot
p_bray <- ggplot(
  ord_df,
  aes(x = Axis.1, y = Axis.2, color = Timepoint)
) +
  geom_point(size = 3, alpha = 0.9) +
  theme_bw(base_size = 13) +
  labs(
    title = "Bray–Curtis beta diversity (abundance-weighted)",
    x = paste0(
      "PCoA1 (",
      round(ord_bray$values$Relative_eig[1] * 100, 1),
      "%)"
    ),
    y = paste0(
      "PCoA2 (",
      round(ord_bray$values$Relative_eig[2] * 100, 1),
      "%)"
    ),
    color = "Timepoint"
  )

p_bray


