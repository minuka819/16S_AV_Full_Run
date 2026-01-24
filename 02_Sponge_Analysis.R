### --- CHAPTER 2 --- ###

### --- SPONGE ONLY ANALYSIS --- ###

PS_16S_Sponge

taxa_df <- get_taxa(PS_16S_Sponge)


PS_sponge <- prune_taxa(taxa_sums(PS_16S_Sponge) > 0, PS_16S_Sponge)


PS_sponge <- subset_taxa(PS_sponge, Order != "Chloroplast" & Family != "Mitochondria")
PS_sponge <- prune_taxa(taxa_sums(PS_sponge) > 0, PS_sponge)

PS_sponge_rare <- subset_samples(PS_rare, Sample_Type == "Sponge")
PS_sponge_rare <- prune_taxa(taxa_sums(PS_sponge_rare) > 0, PS_sponge_rare)


library(phyloseq)
library(ggplot2)
library(patchwork)

# ---- Settings ----
COLORS <- c("Leaf" = "#028A0F", "Sponge" = "#02A3D3")  # kept for consistency (not used here)
TANK_COLORS <- c(
  "Existing_tank" = "#7A7A7A",  # neutral grey (baseline/commercial)
  "Control"       = "#F28E2B",  # orange
  "Inoculated"    = "#4E79A7"   # blue
)
time_levels <- c("Pre_inoculation", "1_week_post_inoculation", "Harvest")

# ---- Use your rarefied sponge object ----
ps <- PS_sponge_rare   # <-- rename if your object is named differently

# Keep only the 3 timepoints (drops samples with NA/other timepoints)
ps <- subset_samples(ps, Timepoint %in% time_levels)
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

# Make factors with explicit ordering
sd <- data.frame(sample_data(ps))
sd$Timepoint <- factor(sd$Timepoint, levels = time_levels)
sd$Tank_Type <- factor(sd$Tank_Type, levels = c("Existing_tank", "Control", "Inoculated"))
sample_data(ps) <- sample_data(sd)

# ---- Alpha diversity table ----
alpha_df <- estimate_richness(ps, measures = c("Observed", "Shannon"))
alpha_df$Sample <- rownames(alpha_df)
alpha_df$Timepoint <- sample_data(ps)$Timepoint
alpha_df$Tank_Type <- sample_data(ps)$Tank_Type

# ---- Plots (Observed + Shannon) ----
p_obs <- ggplot(alpha_df, aes(x = Tank_Type, y = Observed, fill = Tank_Type)) +
  geom_boxplot() +
  geom_jitter(width = 0.15, alpha = 0.7, size = 2) +
  facet_wrap(~ Timepoint, scales = "free_x") +
  scale_fill_manual(values = TANK_COLORS) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(title = "Observed ASVs (rarefied sponge samples)", x = NULL, y = "Observed")

p_sha <- ggplot(alpha_df, aes(x = Tank_Type, y = Shannon, fill = Tank_Type)) +
  geom_boxplot() +
  geom_jitter(width = 0.15, alpha = 0.7, size = 2) +
  facet_wrap(~ Timepoint, scales = "free_x") +
  scale_fill_manual(values = TANK_COLORS) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(title = "Shannon diversity (rarefied sponge samples)", x = NULL, y = "Shannon")

p_obs + p_sha

# ---- Stats: within each timepoint (3-group comparison) ----
# 1) Kruskal-Wallis across Tank_Type, separately per timepoint
kw_results <- by(alpha_df, alpha_df$Timepoint, function(d) {
  data.frame(
    Metric = c("Observed", "Shannon"),
    P_value = c(
      kruskal.test(Observed ~ Tank_Type, data = d)$p.value,
      kruskal.test(Shannon  ~ Tank_Type, data = d)$p.value
    )
  )
})
kw_results

# 2) Pairwise Wilcoxon (if you want post-hoc comparisons), per timepoint
pairwise_results <- by(alpha_df, alpha_df$Timepoint, function(d) {
  list(
    Observed = pairwise.wilcox.test(d$Observed, d$Tank_Type, p.adjust.method = "BH"),
    Shannon  = pairwise.wilcox.test(d$Shannon,  d$Tank_Type, p.adjust.method = "BH")
  )
})
pairwise_results

bray_dist <- phyloseq::distance(PS_sponge, method = "bray")

ord_bray <- ordinate(PS_sponge, method = "PCoA", distance = bray_dist)

TANK_COLORS <- c(
  "Existing_tank" = "#7A7A7A",
  "Control"       = "#F28E2B",
  "Inoculated"    = "#4E79A7"
)


library(ggplot2)

p_bray <- plot_ordination(
  PS_sponge,
  ord_bray,
  color = "Tank_Type",
  shape = "Timepoint"
) +
  geom_point(size = 3, alpha = 0.85) +
  scale_color_manual(values = TANK_COLORS) +
  theme_bw() +
  labs(
    title = "Bray–Curtis PCoA (sponge samples)",
    color = "Tank Type",
    shape = "Timepoint"
  )

p_bray


library(vegan)

# Extract metadata
meta <- data.frame(sample_data(PS_sponge))

# PERMANOVA: main effects + interaction
adonis_res <- adonis2(
  bray_dist ~ Tank_Type * Timepoint,
  data = meta,
  permutations = 999
)

adonis_res

plot_ordination(
  PS_sponge,
  ord_bray,
  color = "Tank_Type"
) +
  geom_point(size = 3, alpha = 0.85) +
  scale_color_manual(values = TANK_COLORS) +
  theme_bw() +
  labs(
    title = "Bray–Curtis PCoA by tank type (sponge samples)",
    color = "Tank Type"
  )

plot_ordination(
  PS_sponge,
  ord_bray,
  color = "Timepoint"
) +
  geom_point(size = 3, alpha = 0.85) +
  theme_bw() +
  labs(
    title = "Bray–Curtis PCoA by timepoint (sponge samples)",
    color = "Timepoint"
  )



rank_names(PS_sponge)



library(phyloseq)

# Normalize genus names to character
tax <- as.data.frame(tax_table(PS_sponge))
tax$Genus <- as.character(tax$Genus)

# Subset taxa that match genera of interest
keep_taxa <- taxa_names(PS_sponge)[
  tax$Genus %in% genera_interest
]

PS_sponge_genera <- prune_taxa(keep_taxa, PS_sponge)
PS_sponge_genera <- prune_taxa(taxa_sums(PS_sponge_genera) > 0, PS_sponge_genera)

PS_sponge_genera

PS_actino <- subset_taxa(PS_sponge, Phylum == "Actinobacteriota")


library(dplyr)
library(ggplot2)

# Convert to relative abundance
PS_rel <- transform_sample_counts(PS_sponge_genera, function(x) x / sum(x))

df_rel <- psmelt(PS_rel)

# Clean genus names
df_rel$Genus <- as.character(df_rel$Genus)
df_rel$Genus[df_rel$Genus == "" | is.na(df_rel$Genus)] <- "Unassigned"

# Keep only Control vs Inoculated
df_rel <- df_rel %>%
  filter(Tank_Type %in% c("Control", "Inoculated"))

df_rel$Timepoint <- factor(
  df_rel$Timepoint,
  levels = c(
    "Pre_inoculation",
    "1_week_post_inoculation",
    "Harvest"
  )
)

df_sum <- df_rel %>%
  group_by(Timepoint, Tank_Type, Genus) %>%
  summarise(
    Mean_RA = mean(Abundance),
    .groups = "drop"
  )

TANK_COLORS <- c(
  "Control"    = "#F28E2B",
  "Inoculated" = "#4E79A7"
)

ggplot(df_sum, aes(x = Genus, y = Mean_RA, fill = Tank_Type)) +
  geom_col(position = "dodge") +
  facet_wrap(~ Timepoint, scales = "free_y") +
  scale_fill_manual(values = TANK_COLORS) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    title = "Relative abundance of genera of interest (sponge samples)",
    x = "Genus",
    y = "Mean relative abundance",
    fill = "Tank Type"
  )


