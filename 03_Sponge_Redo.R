
#Rarefying sponge only 
PS_sponge <- subset_samples(PS_16S, Sample_Type == "Sponge")


# Get sequencing depth per sample
sponge_depths <- sample_sums(PS_sponge)

# Inspect depths
summary(sponge_depths)

sponge_depths
view.as.data.frame(sponge_depths)

# Minimum depth
min_depth <- min(sponge_depths)
min_depth

depth_df <- data.frame(
  SampleID = names(sample_sums(PS_sponge)),
  Depth = sample_sums(PS_sponge)
)

# Join metadata (optional but useful for coloring/faceting later)
depth_df <- cbind(
  depth_df,
  sample_data(PS_sponge)[depth_df$SampleID, ]
)


ggplot(depth_df, aes(x = reorder(SampleID, Depth), y = Depth)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  theme_bw(base_size = 13) +
  labs(
    title = "Read depth per sponge sample",
    x = "Sample",
    y = "Sequencing depth (reads)"
  )


#Rarefy

set.seed(123)

PS_sponge_rarefied <- rarefy_even_depth(
  PS_sponge,
  sample.size = 39000,
  rngseed = 123,
  replace = FALSE,
  verbose = TRUE
)

sample_sums(PS_sponge_rarefied)

alpha_df <- estimate_richness(
  PS_sponge_rarefied,
  measures = c("Observed", "Shannon")
)

alpha_df$SampleID <- rownames(alpha_df)

alpha_df <- cbind(
  alpha_df,
  sample_data(PS_sponge_rarefied)[alpha_df$SampleID, ]
)

# Factor ordering (CRITICAL for clean panels)
alpha_df$Tank_Type <- factor(
  alpha_df$Tank_Type,
  levels = names(TANK_COLORS)
)

alpha_df$Timepoint <- factor(
  alpha_df$Timepoint,
  levels = time_levels
)

p_obs <- ggplot(
  alpha_df,
  aes(x = Tank_Type, y = Observed, fill = Tank_Type)
) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.7, size = 2) +
  facet_wrap(~ Timepoint, scales = "free_x") +
  scale_fill_manual(values = TANK_COLORS) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold")
  ) +
  labs(
    title = "Observed ASVs (rarefied sponge samples)",
    x = NULL,
    y = "Observed"
  )

p_obs


p_sha <- ggplot(
  alpha_df,
  aes(x = Tank_Type, y = Shannon, fill = Tank_Type)
) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.7, size = 2) +
  facet_wrap(~ Timepoint, scales = "free_x") +
  scale_fill_manual(values = TANK_COLORS) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold")
  ) +
  labs(
    title = "Shannon diversity (rarefied sponge samples)",
    x = NULL,
    y = "Shannon"
  )

p_sha


p_obs + p_sha


alpha_stats_df <- alpha_df %>%
  dplyr::filter(Tank_Type %in% c("Control", "Inoculated"))

run_wilcox_by_time <- function(df, metric) {
  df %>%
    dplyr::group_by(Timepoint) %>%
    dplyr::summarise(
      metric = metric,
      p_value = wilcox.test(
        formula = as.formula(paste(metric, "~ Tank_Type")),
        data = cur_data()
      )$p.value
    )
}


obs_stats <- run_wilcox_by_time(alpha_stats_df, "Observed")
obs_stats

sha_stats <- run_wilcox_by_time(alpha_stats_df, "Shannon")
sha_stats

obs_stats$p_adj <- p.adjust(obs_stats$p_value, method = "BH")
sha_stats$p_adj <- p.adjust(sha_stats$p_value, method = "BH")

obs_stats
sha_stats


library(vegan)

# Distance matrix
bray_dist <- phyloseq::distance(PS_sponge, method = "bray")

ord_bray <- ordinate(
  PS_sponge,
  method = "PCoA",
  distance = bray_dist
)

p_bray <- plot_ordination(
  PS_sponge,
  ord_bray,
  color = "Tank_Type",
  shape = "Timepoint"
) +
  geom_point(size = 3, alpha = 0.85) +
  scale_color_manual(values = TANK_COLORS) +
  theme_bw(base_size = 13) +
  labs(
    title = "Bray–Curtis PCoA (sponge samples)",
    color = "Tank Type",
    shape = "Timepoint"
  )

p_bray


library(vegan)

# Metadata
meta_bray <- data.frame(sample_data(PS_sponge))

# Bray–Curtis distance
bray_dist <- phyloseq::distance(PS_sponge, method = "bray")

# PERMANOVA
adonis_bray <- adonis2(
  bray_dist ~ Tank_Type + Timepoint,
  data = meta_bray,
  permutations = 999
)

adonis_bray

disp_bray <- betadisper(bray_dist, meta_bray$Tank_Type)
anova(disp_bray)


# Jaccard distance (presence/absence)
jacc_dist <- phyloseq::distance(
  PS_sponge_rarefied,
  method = "jaccard",
  binary = TRUE
)


ord_jacc <- ordinate(
  PS_sponge_rarefied,
  method = "PCoA",
  distance = jacc_dist
)

p_jacc <- plot_ordination(
  PS_sponge_rarefied,
  ord_jacc,
  color = "Tank_Type",
  shape = "Timepoint"
) +
  geom_point(size = 3, alpha = 0.85) +
  scale_color_manual(values = TANK_COLORS) +
  theme_bw(base_size = 13) +
  labs(
    title = "Jaccard PCoA (presence/absence, rarefied sponge samples)",
    color = "Tank Type",
    shape = "Timepoint"
  )

p_jacc


meta_jacc <- data.frame(sample_data(PS_sponge_rarefied))

adonis_jacc <- adonis2(
  jacc_dist ~ Tank_Type + Timepoint,
  data = meta_jacc,
  permutations = 999
)

adonis_jacc


disp_jacc <- betadisper(jacc_dist, meta_jacc$Tank_Type)
anova(disp_jacc)


#Actinobacteria phylum


PS_actino <- subset_taxa(PS_sponge_rare, Phylum == "Actinobacteriota")

PS_actino_rel <- transform_sample_counts(
  PS_actino,
  function(x) x / sum(x)
)


PS_actino_phylum <- tax_glom(
  PS_actino_rel,
  taxrank = "Phylum"
)


actino_df <- psmelt(PS_actino_phylum)


actino_df$Tank_Type <- factor(
  actino_df$Tank_Type,
  levels = names(TANK_COLORS)
)

actino_df$Timepoint <- factor(
  actino_df$Timepoint,
  levels = time_levels
)


actino_summary <- actino_df %>%
  dplyr::group_by(Tank_Type, Timepoint) %>%
  dplyr::summarise(
    mean_abund = mean(Abundance),
    se_abund = sd(Abundance) / sqrt(n()),
    .groups = "drop"
  )

ggplot(
  actino_summary,
  aes(x = Tank_Type, y = mean_abund, fill = Tank_Type)
) +
  geom_col(position = position_dodge(width = 0.7)) +
  geom_errorbar(
    aes(ymin = mean_abund - se_abund, ymax = mean_abund + se_abund),
    width = 0.2,
    position = position_dodge(width = 0.7)
  ) +
  facet_wrap(~ Timepoint) +
  scale_fill_manual(values = TANK_COLORS) +
  theme_bw(base_size = 13) +
  theme(legend.position = "none") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title = "Mean relative abundance of Actinobacteria",
    x = NULL,
    y = "Relative abundance (%)"
  )


PS_actino <- subset_taxa(
  PS_sponge,
  Phylum == "Actinobacteriota" | Phylum == "Actinobacteria"
)


PS_actino_rel <- transform_sample_counts(
  PS_actino,
  function(x) x / sum(x)
)


PS_actino_genus <- tax_glom(
  PS_actino_rel,
  taxrank = "Genus",
  NArm = TRUE
)

actino_genus_df <- psmelt(PS_actino_genus)


actino_genus_df$Tank_Type <- factor(
  actino_genus_df$Tank_Type,
  levels = names(TANK_COLORS)
)

actino_genus_df$Timepoint <- factor(
  actino_genus_df$Timepoint,
  levels = time_levels
)


top_genera <- actino_genus_df %>%
  dplyr::group_by(Genus) %>%
  dplyr::summarise(total_abundance = sum(Abundance)) %>%
  dplyr::arrange(desc(total_abundance)) %>%
  dplyr::slice_head(n = 10) %>%
  dplyr::pull(Genus)


actino_genus_df$Genus_plot <- ifelse(
  actino_genus_df$Genus %in% top_genera,
  actino_genus_df$Genus,
  "Other"
)

actino_genus_summary <- actino_genus_df %>%
  dplyr::group_by(Timepoint, Tank_Type, Genus_plot) %>%
  dplyr::summarise(
    mean_abund = mean(Abundance),
    .groups = "drop"
  )


p_actino_genus <- ggplot(
  actino_genus_summary,
  aes(x = Tank_Type, y = mean_abund, fill = Genus_plot)
) +
  geom_col() +
  facet_wrap(~ Timepoint) +
  theme_bw(base_size = 13) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title = "Actinobacteria genus-level composition (sponge samples)",
    x = NULL,
    y = "Mean relative abundance (%)",
    fill = "Genus"
  )

p_actino_genus


### NON STACKED


actino_genus_summary$Genus_plot <- factor(
  actino_genus_summary$Genus_plot,
  levels = c(sort(unique(actino_genus_summary$Genus_plot[actino_genus_summary$Genus_plot != "Other"])), "Other")
)



p_actino_genus_grouped <- ggplot(
  actino_genus_summary,
  aes(x = Tank_Type, y = mean_abund, fill = Tank_Type)
) +
  geom_col(position = position_dodge(width = 0.7)) +
  facet_grid(Genus_plot ~ Timepoint, scales = "free_y") +
  scale_fill_manual(values = TANK_COLORS) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    strip.text.y = element_text(angle = 0)
  ) +
  labs(
    title = "Actinobacteria genus-level relative abundance (sponge samples)",
    x = NULL,
    y = "Mean relative abundance (%)"
  )

p_actino_genus_grouped


###RETRY

PS_actino <- subset_taxa(
  PS_sponge,
  Phylum %in% c("Actinobacteria", "Actinobacteriota")
)


PS_actino_rel <- transform_sample_counts(
  PS_actino,
  function(x) x / sum(x)
)

PS_actino_genus <- tax_glom(
  PS_actino_rel,
  taxrank = "Genus",
  NArm = TRUE
)


actino_df <- psmelt(PS_actino_genus)


actino_df <- actino_df %>%
  dplyr::filter(
    Tank_Type %in% c("Control", "Inoculated"),
    !is.na(Genus),
    Genus != ""
  )

actino_df$Timepoint <- factor(
  actino_df$Timepoint,
  levels = time_levels
)

actino_df$Tank_Type <- factor(
  actino_df$Tank_Type,
  levels = c("Control", "Inoculated")
)



actino_summary <- actino_df %>%
  dplyr::group_by(Timepoint, Tank_Type, Genus) %>%
  dplyr::summarise(
    mean_abund = mean(Abundance),
    .groups = "drop"
  )

actino_summary <- actino_summary %>%
  dplyr::group_by(Genus) %>%
  dplyr::filter(max(mean_abund) > 0.001) %>%  # 0.1%
  ungroup()


p_actino_genus <- ggplot(
  actino_summary,
  aes(x = Genus, y = mean_abund, fill = Tank_Type)
) +
  geom_col(position = position_dodge(width = 0.7)) +
  facet_wrap(~ Timepoint) +
  scale_fill_manual(values = TANK_COLORS) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  ) +
  labs(
    title = "Relative abundance of Actinobacteria genera (sponge samples)",
    x = "Genus",
    y = "Mean relative abundance (%)",
    fill = "Tank Type"
  )

p_actino_genus


meta_all<-data.frame(sample_data(PS_16S))


library(phyloseq)
library(dplyr)
library(ggplot2)
library(scales)

# -------------------------------
# Actinobacteria genus-level plot
# -------------------------------

# 1. Subset sponge samples to Actinobacteria only
PS_actino <- subset_taxa(
  PS_sponge,
  Phylum %in% c("Actinobacteria", "Actinobacteriota")
)

# 2. Convert to relative abundance
PS_actino_rel <- transform_sample_counts(
  PS_actino,
  function(x) x / sum(x)
)

# 3. Agglomerate at genus level
PS_actino_genus <- tax_glom(
  PS_actino_rel,
  taxrank = "Genus",
  NArm = TRUE
)

# 4. Melt to dataframe
actino_df <- psmelt(PS_actino_genus)

# 5. Clean taxonomy + metadata
actino_df <- actino_df %>%
  filter(
    !is.na(Genus),
    Genus != ""
  )

actino_df$Timepoint <- factor(
  actino_df$Timepoint,
  levels = c("Pre_inoculation", "1_week_post_inoculation", "Harvest")
)

# 6. Mean relative abundance per genus × timepoint
actino_summary <- actino_df %>%
  group_by(Genus, Timepoint) %>%
  summarise(
    mean_abund = mean(Abundance),
    .groups = "drop"
  )

# 7. Select top 20 Actinobacteria genera (overall)
top20_genera <- actino_summary %>%
  group_by(Genus) %>%
  summarise(total_abund = sum(mean_abund)) %>%
  arrange(desc(total_abund)) %>%
  slice_head(n = 20) %>%
  pull(Genus)

actino_summary_top20 <- actino_summary %>%
  filter(Genus %in% top20_genera)

# Preserve genus order
actino_summary_top20$Genus <- factor(
  actino_summary_top20$Genus,
  levels = top20_genera
)

# 8. Timepoint colors (matches your example)
TIME_COLORS <- c(
  "Pre_inoculation" = "#E15759",
  "1_week_post_inoculation" = "#2CA02C",
  "Harvest" = "#4E79A7"
)

# 9. Plot
ggplot(
  actino_summary_top20,
  aes(x = Genus, y = mean_abund, fill = Timepoint)
) +
  geom_col(position = position_dodge(width = 0.7)) +
  scale_fill_manual(values = TIME_COLORS) +
  scale_y_continuous() +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  ) +
  labs(
    title = "Top 20 Actinobacteria genera (relative abundance, sponge samples)",
    x = "Genus",
    y = "Mean relative abundance (%)",
    fill = "Timepoint"
  )

library(phyloseq)
library(dplyr)
library(ggplot2)
library(scales)

# -------------------------------
# Actinobacteria genus-level plot
# Control vs Inoculated + Time
# -------------------------------

# 1. Subset Actinobacteria (sponge only assumed)
PS_actino <- subset_taxa(
  PS_sponge,
  Phylum %in% c("Actinobacteria", "Actinobacteriota")
)

# 2. Relative abundance
PS_actino_rel <- transform_sample_counts(
  PS_actino,
  function(x) x / sum(x)
)

# 3. Agglomerate at genus
PS_actino_genus <- tax_glom(
  PS_actino_rel,
  taxrank = "Genus",
  NArm = TRUE
)

# 4. Melt
actino_df <- psmelt(PS_actino_genus)

# 5. Clean metadata
actino_df <- actino_df %>%
  filter(
    Tank_Type %in% c("Control", "Inoculated"),
    !is.na(Genus),
    Genus != ""
  )

actino_df$Timepoint <- factor(
  actino_df$Timepoint,
  levels = c("Pre_inoculation", "1_week_post_inoculation", "Harvest")
)

actino_df$Tank_Type <- factor(
  actino_df$Tank_Type,
  levels = c("Control", "Inoculated")
)

# 6. Mean abundance per Genus × Tank × Time
actino_summary <- actino_df %>%
  group_by(Genus, Tank_Type, Timepoint) %>%
  summarise(
    mean_abund = mean(Abundance),
    .groups = "drop"
  )

# 7. Top 20 genera overall
top20_genera <- actino_summary %>%
  group_by(Genus) %>%
  summarise(total_abund = sum(mean_abund)) %>%
  arrange(desc(total_abund)) %>%
  slice_head(n = 20) %>%
  pull(Genus)

actino_summary <- actino_summary %>%
  filter(Genus %in% top20_genera)

actino_summary$Genus <- factor(
  actino_summary$Genus,
  levels = top20_genera
)

# 8. Plot
ggplot(
  actino_summary,
  aes(x = Genus, y = mean_abund, fill = Tank_Type)
) +
  geom_col(position = position_dodge(width = 0.7)) +
  facet_wrap(~ Timepoint) +
  scale_fill_manual(values = c(
    "Control" = "#F28E2B",
    "Inoculated" = "#4E79A7"
  )) +
  scale_y_continuous() +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  ) +
  labs(
    title = "Top 20 Actinobacteria genera (relative abundance, sponge samples)",
    x = "Genus",
    y = "Mean relative abundance (%)",
    fill = "Tank Type"
  )


library(phyloseq)
library(dplyr)
library(ggplot2)
library(scales)

# -------------------------------
# Actinobacteria: Top 20 genera
# (bars = mean relative abundance)
# (x labels = total read counts)
# -------------------------------
library(phyloseq)
library(dplyr)
library(ggplot2)

# -------------------------------
# Actinobacteria genus plot:
# fill = Tank_Type, facet = Timepoint
# + read counts under each Genus label
# -------------------------------

# 1) Subset Actinobacteria (sponge samples assumed in PS_sponge)
PS_actino <- subset_taxa(
  PS_sponge,
  Phylum %in% c("Actinobacteria", "Actinobacteriota")
)

# 2) Genus-level RAW counts (for label counts)
PS_actino_genus_raw <- tax_glom(PS_actino, taxrank = "Genus", NArm = TRUE)

# Total reads per genus across ALL sponge samples
genus_counts <- data.frame(
  Genus = as.character(tax_table(PS_actino_genus_raw)[, "Genus"]),
  total_reads = taxa_sums(PS_actino_genus_raw)
) %>%
  filter(!is.na(Genus), Genus != "") %>%
  group_by(Genus) %>%
  summarise(total_reads = sum(total_reads), .groups = "drop")

# 3) Genus-level RELATIVE abundance (for plotting)
PS_actino_rel <- transform_sample_counts(PS_actino, function(x) x / sum(x))
PS_actino_genus_rel <- tax_glom(PS_actino_rel, taxrank = "Genus", NArm = TRUE)

actino_df <- psmelt(PS_actino_genus_rel) %>%
  filter(
    Tank_Type %in% c("Control", "Inoculated"),
    !is.na(Genus),
    Genus != ""
  )

# Factor ordering (match your pipeline)
actino_df$Timepoint <- factor(
  actino_df$Timepoint,
  levels = c("Pre_inoculation", "1_week_post_inoculation", "Harvest")
)
actino_df$Tank_Type <- factor(actino_df$Tank_Type, levels = c("Control", "Inoculated"))

# 4) Mean relative abundance per Genus × Tank_Type × Timepoint
actino_summary <- actino_df %>%
  group_by(Genus, Tank_Type, Timepoint) %>%
  summarise(mean_abund = mean(Abundance), .groups = "drop")

# 5) Top 20 genera overall (based on total mean abundance across groups/time)
top20_genera <- actino_summary %>%
  group_by(Genus) %>%
  summarise(total_abund = sum(mean_abund), .groups = "drop") %>%
  arrange(desc(total_abund)) %>%
  slice_head(n = 20) %>%
  pull(Genus)

actino_summary_top20 <- actino_summary %>%
  filter(Genus %in% top20_genera) %>%
  left_join(genus_counts, by = "Genus")

# Keep genus order consistent
actino_summary_top20$Genus <- factor(actino_summary_top20$Genus, levels = top20_genera)

# 6) Build x-axis labels with read counts
# (use base R formatting to avoid scales dependency issues)
label_map <- setNames(
  paste0(
    top20_genera,
    "\n(n=", format(genus_counts$total_reads[match(top20_genera, genus_counts$Genus)],
                    big.mark = ",", scientific = FALSE),
    ")"
  ),
  top20_genera
)

# ---- Make a label column (mean %; optional add reads) ----
actino_summary_top20 <- actino_summary_top20 %>%
  mutate(
    pct_label = sprintf("%.1f%%", mean_abund * 100),
    reads_label = format(total_reads, big.mark = ",", scientific = FALSE),
    bar_label = paste0(pct_label, "\n(n=", reads_label, ")")  # or just pct_label
  )

library(dplyr)
library(ggplot2)

lab_thresh <- 0.01   # 1% relative abundance

plot_df <- actino_summary_top20 %>%
  dplyr::mutate(
    read_label = ifelse(
      mean_abund >= lab_thresh,
      format(total_reads, big.mark = ",", scientific = FALSE),
      ""
    )
  )

ggplot(
  plot_df,
  aes(x = Genus, y = mean_abund, fill = Tank_Type)
) +
  geom_col(position = position_dodge(width = 0.7)) +
  geom_text(
    aes(label = read_label),
    position = position_dodge(width = 0.7),
    angle = 90,          # <-- vertical text
    vjust = -0.1,
    size = 3
  ) +
  facet_wrap(~ Timepoint) +
  scale_fill_manual(values = c(
    "Control" = "#F28E2B",
    "Inoculated" = "#4E79A7"
  )) +
  scale_y_continuous(
    labels = function(x) paste0(x * 100, "%"),
    expand = expansion(mult = c(0.02, 0.25))  # headroom for vertical labels
  ) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  ) +
  labs(
    title = "Top 20 Actinobacteria genera (relative abundance; read counts shown)",
    x = "Genus",
    y = "Mean relative abundance (%)",
    fill = "Tank Type"
  )


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


library(phyloseq)
library(dplyr)
library(ggplot2)

# ---------------------------------------
# Genera of interest – relative abundance
# (same plot logic as Actinobacteria)
# ---------------------------------------

PS_sponge_no_control <- subset_samples(
  PS_sponge,
  Tank_Type != "N/A"
)


PS_sponge_control <- subset_samples(
  PS_sponge,
  Sample_Number == "73"
)

PS_sponge_control

# 1) Subset genera of interest (sponge samples assumed)
PS_goi <- subset_taxa(
  PS_sponge_no_control,
  Genus %in% genera_of_interest
)

# 2) RAW counts at genus level (for read labels)
PS_goi_genus_raw <- tax_glom(PS_goi, taxrank = "Genus", NArm = TRUE)

genus_counts <- data.frame(
  Genus = as.character(tax_table(PS_goi_genus_raw)[, "Genus"]),
  total_reads = taxa_sums(PS_goi_genus_raw)
) %>%
  filter(!is.na(Genus), Genus != "") %>%
  group_by(Genus) %>%
  summarise(total_reads = sum(total_reads), .groups = "drop")

# 3) Relative abundance (for plotting)
PS_goi_rel <- transform_sample_counts(PS_goi, function(x) x / sum(x))
PS_goi_genus_rel <- tax_glom(PS_goi_rel, taxrank = "Genus", NArm = TRUE)

goi_df <- psmelt(PS_goi_genus_rel) %>%
  filter(
    Tank_Type %in% c("Control", "Inoculated"),
    !is.na(Genus),
    Genus != ""
  )

# Factor ordering
goi_df$Timepoint <- factor(
  goi_df$Timepoint,
  levels = c("Pre_inoculation", "1_week_post_inoculation", "Harvest")
)

goi_df$Tank_Type <- factor(
  goi_df$Tank_Type,
  levels = c("Control", "Inoculated")
)

goi_df$Genus <- factor(goi_df$Genus, levels = genera_of_interest)

# 4) Mean relative abundance per Genus × Tank × Time
goi_summary <- goi_df %>%
  group_by(Genus, Tank_Type, Timepoint) %>%
  summarise(mean_abund = mean(Abundance), .groups = "drop") %>%
  left_join(genus_counts, by = "Genus")

# 5) Add vertical read-count labels (thresholded)
lab_thresh <- 0.01  # 1%

goi_summary <- goi_summary %>%
  mutate(
    read_label = ifelse(
      mean_abund >= lab_thresh,
      format(total_reads, big.mark = ",", scientific = FALSE),
      ""
    )
  )

# 6) Plot (EXACT same style as Actinobacteria)
ggplot(
  goi_summary,
  aes(x = Genus, y = mean_abund, fill = Tank_Type)
) +
  geom_col(position = position_dodge(width = 0.7)) +
  geom_text(
    aes(label = read_label),
    position = position_dodge(width = 0.7),
    angle = 90,
    vjust = -0.1,
    size = 3
  ) +
  facet_wrap(~ Timepoint) +
  scale_fill_manual(values = c(
    "Control" = "#F28E2B",
    "Inoculated" = "#4E79A7"
  )) +
  scale_y_continuous(
    labels = function(x) paste0(x * 100, "%"),
    expand = expansion(mult = c(0.02, 0.25))
  ) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  ) +
  labs(
    title = "Relative abundance of genera of interest (sponge samples)",
    x = "Genus",
    y = "Mean relative abundance (%)",
    fill = "Tank Type"
  )

meta_sponge <- sample_data(PS_sponge)


library(phyloseq)
library(ggplot2)

# Convert to relative abundance
PS_ctrl_rel <- transform_sample_counts(
  PS_sponge_control,
  function(x) x / sum(x)
)

# Collapse at Phylum level
PS_ctrl_phylum <- tax_glom(
  PS_ctrl_rel,
  taxrank = "Phylum"
)

# Quick stacked bar plot
plot_bar(
  PS_ctrl_phylum,
  fill = "Phylum"
) +
  theme_bw(base_size = 13) +
  labs(
    title = "Phylum-level relative abundance (PS_sponge_control)",
    y = "Relative abundance"
  )

# Agglomerate at Genus
PS_ctrl_genus <- tax_glom(
  PS_ctrl_rel,
  taxrank = "Genus",
  NArm = TRUE
)

# Keep top 10 genera
top10_genera <- names(sort(taxa_sums(PS_ctrl_genus), decreasing = TRUE))[1:10]

PS_ctrl_genus_top10 <- prune_taxa(
  top10_genera,
  PS_ctrl_genus
)

# Plot
plot_bar(
  PS_ctrl_genus_top10,
  fill = "Genus"
) +
  theme_bw(base_size = 13) +
  labs(
    title = "Top 10 genera (relative abundance) – PS_sponge_control",
    y = "Relative abundance"
  )
