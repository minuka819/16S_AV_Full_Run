### --- CHAPTER 2 --- ###

### --- SPONGE ONLY ANALYSIS --- ###

PS_16S_Sponge

taxa_df <- as.data.frame(tax_table(PS_sponge))


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


p_obs


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

bray_dist <- phyloseq::distance(PS_sponge_rare, method = "bray")

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
meta <- data.frame(sample_data(PS_sponge_rare))

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


#### - JACCARD BETA - ####


library(phyloseq)
library(vegan)
library(ggplot2)

# ---- 1) Subset rarefied object to sponge + experimental tanks ----
ps <- subset_samples(PS_rare, Sample_Type == "Sponge" & Tank_Type %in% c("Control", "Inoculated"))
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

# ---- 2) Build presence/absence OTU table safely (preserve taxa names) ----
otu <- as(otu_table(ps), "matrix")

# Ensure rows = taxa
if (!taxa_are_rows(ps)) otu <- t(otu)

otu_pa <- (otu > 0) * 1
rownames(otu_pa) <- rownames(otu)
colnames(otu_pa) <- colnames(otu)

otu_pa <- otu_table(otu_pa, taxa_are_rows = TRUE)

# ---- 3) Rebuild phyloseq object (drop tree to avoid taxa mismatch issues) ----
tax <- tax_table(ps)
sam <- sample_data(ps)

PS_pa <- phyloseq(otu_pa, tax, sam)
PS_pa <- prune_taxa(taxa_sums(PS_pa) > 0, PS_pa)  # just to be safe

# ---- 4) Jaccard distance + PCoA ----
jac_dist <- phyloseq::distance(PS_pa, method = "jaccard", binary = TRUE)
ord_jac  <- ordinate(PS_pa, method = "PCoA", distance = jac_dist)

TANK_COLORS <- c("Control" = "#F28E2B", "Inoculated" = "#4E79A7")

p_jac <- plot_ordination(PS_pa, ord_jac, color = "Tank_Type", shape = "Timepoint") +
  geom_point(size = 3, alpha = 0.85) +
  scale_color_manual(values = TANK_COLORS) +
  theme_bw() +
  labs(
    title = "Jaccard PCoA (presence/absence) - rarefied sponge samples",
    color = "Tank Type",
    shape = "Timepoint"
  )

p_jac

# ---- 5) PERMANOVA (adonis2) ----
meta <- data.frame(sample_data(PS_pa))
meta$Timepoint <- factor(meta$Timepoint,
                         levels = c("Pre_inoculation", "1_week_post_inoculation", "Harvest"))

adonis_overall <- adonis2(jac_dist ~ Tank_Type * Timepoint, data = meta, permutations = 999)
adonis_terms   <- adonis2(jac_dist ~ Tank_Type * Timepoint, data = meta, permutations = 999, by = "terms")

adonis_overall
adonis_terms

# ---- 6) Dispersion check (important) ----
disp_tank <- betadisper(jac_dist, meta$Tank_Type)
anova(disp_tank)





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

PS_actino <- subset_taxa(PS_sponge, Phylum == "Actinobacteriota")
PS_actino <- prune_taxa(taxa_sums(PS_actino) > 0, PS_actino)

# Convert to relative abundance within each sample
PS_actino_rel <- transform_sample_counts(
  PS_actino,
  function(x) x / sum(x)
)


library(dplyr)

df_actino <- psmelt(PS_actino_rel)

# Keep experimental tanks only
df_actino <- df_actino %>%
  filter(Tank_Type %in% c("Control", "Inoculated"))

# Enforce biological time order
df_actino$Timepoint <- factor(
  df_actino$Timepoint,
  levels = c(
    "Pre_inoculation",
    "1_week_post_inoculation",
    "Harvest"
  )
)

df_actino_sum <- df_actino %>%
  group_by(Timepoint, Tank_Type) %>%
  summarise(
    Mean_RA = mean(Abundance),
    .groups = "drop"
  )

library(ggplot2)

TANK_COLORS <- c(
  "Control"    = "#F28E2B",
  "Inoculated" = "#4E79A7"
)

ggplot(df_actino_sum,
       aes(x = Timepoint, y = Mean_RA, fill = Tank_Type)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = TANK_COLORS) +
  theme_bw() +
  labs(
    title = "Actinobacteria relative abundance in sponge samples",
    x = "Timepoint",
    y = "Mean relative abundance",
    fill = "Tank Type"
  )


library(phyloseq)
library(dplyr)
library(ggplot2)

# 1) Build Actinobacteria-only object (if not already created)
PS_actino <- subset_taxa(PS_sponge_rare, Phylum == "Actinobacteriota")
PS_actino <- prune_taxa(taxa_sums(PS_actino) > 0, PS_actino)

# 2) Absolute Actinobacteria reads per sample (counts)
actino_counts <- data.frame(
  Sample = sample_names(PS_actino),
  Actino_Reads = sample_sums(PS_actino),
  Tank_Type = sample_data(PS_actino)$Tank_Type,
  Timepoint = sample_data(PS_actino)$Timepoint
)

# Keep experimental tanks only
actino_counts <- actino_counts %>%
  filter(Tank_Type %in% c("Control", "Inoculated"))

# Order timepoints
actino_counts$Timepoint <- factor(
  actino_counts$Timepoint,
  levels = c("Pre_inoculation", "1_week_post_inoculation", "Harvest")
)

# Colors
TANK_COLORS <- c("Control" = "#F28E2B", "Inoculated" = "#4E79A7")

# 3) Plot: mean ± SE bars + individual points
ggplot(actino_counts, aes(x = Timepoint, y = Actino_Reads, fill = Tank_Type)) +
  stat_summary(fun = mean, geom = "bar", position = position_dodge(width = 0.75), width = 0.65) +
  stat_summary(fun.data = mean_se, geom = "errorbar",
               position = position_dodge(width = 0.75), width = 0.2) +
  geom_point(aes(color = Tank_Type),
             position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.75),
             size = 2, alpha = 0.85) +
  scale_fill_manual(values = TANK_COLORS) +
  scale_color_manual(values = TANK_COLORS) +
  theme_bw() +
  labs(
    title = "Actinobacteria absolute abundance (read counts) in sponge samples",
    x = "Timepoint",
    y = "Actinobacteria read counts",
    fill = "Tank Type",
    color = "Tank Type"
  )

PS_actino_g <- tax_glom(PS_actino, taxrank = "Genus", NArm = FALSE)


PS_actino_rel <- transform_sample_counts(
  PS_actino_g,
  function(x) x / sum(x)
)

top20_genera <- names(
  sort(taxa_sums(PS_actino_rel), decreasing = TRUE)
)[1:20]

PS_actino_top20 <- prune_taxa(top20_genera, PS_actino_rel)


library(dplyr)
library(ggplot2)

df_actino <- psmelt(PS_actino_top20)

# Clean genus names
df_actino$Genus <- as.character(df_actino$Genus)
df_actino$Genus[df_actino$Genus == "" | is.na(df_actino$Genus)] <- "Unassigned"

# Enforce biological time order
df_actino$Timepoint <- factor(
  df_actino$Timepoint,
  levels = c(
    "Pre_inoculation",
    "1_week_post_inoculation",
    "Harvest"
  )
)

df_actino_sum <- df_actino %>%
  group_by(Timepoint, Genus) %>%
  summarise(
    Mean_RA = mean(Abundance),
    .groups = "drop"
  )

ggplot(df_actino_sum, aes(x = Genus, y = Mean_RA, fill = Timepoint)) +
  geom_col(position = "dodge") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    title = "Top 20 Actinobacteria genera (relative abundance, sponge samples)",
    x = "Genus",
    y = "Mean relative abundance",
    fill = "Timepoint"
  )



