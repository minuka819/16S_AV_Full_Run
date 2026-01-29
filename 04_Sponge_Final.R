###############################################################################
# Relative abundance plot by Sample_Type (Top 30 taxa, no collapsing)
#
# Purpose:
#   - Visualize relative abundance of the top 30 most abundant taxa
#     aggregated by Sample_Type
#   - Uses NON-rarefied data
#   - No collapsing into "Other"; only taxa with non-zero abundance are shown
#
# Input:
#   - PS_16S (phyloseq object with all samples)
#
# Notes:
#   - Relative abundance plots are descriptive
#   - Fewer than 30 taxa may appear if others have zero abundance
###############################################################################

library(phyloseq)
library(tidyverse)
library(viridis)
library(ggpubr) #stat
library(Polychrome)


# 1) Transform to relative abundance
PS_rel <- transform_sample_counts(
  PS_16S,
  function(x) x / sum(x)
)

# 2) Identify top 30 taxa globally
top30_taxa <- names(sort(taxa_sums(PS_rel), decreasing = TRUE))[1:50]

# 3) Prune phyloseq object to top 30 taxa only
PS_rel_top30 <- prune_taxa(top30_taxa, PS_rel)

# 4) Melt to long format
rel_df <- psmelt(PS_rel_top30)

# 5) Aggregate within Sample_Type
rel_df_sum <- rel_df %>%
  group_by(Sample_Type, Genus) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

# 6) Re-normalize so each Sample_Type sums to 1
rel_df_sum <- rel_df_sum %>%
  group_by(Sample_Type) %>%
  mutate(Abundance = Abundance / sum(Abundance)) %>%
  ungroup()

n_taxa <- length(unique(rel_df_sum$Genus))

distinct_colors <- palette36.colors(n_taxa)
names(distinct_colors) <- levels(rel_df_sum$Genus)


# 7) Plot
ggplot(rel_df_sum, aes(x = Sample_Type, y = Abundance, fill = Genus)) +
  geom_bar(
    stat = "identity",
    color = "black",
    linewidth = 0.15
  ) +
  scale_fill_manual(values = distinct_colors) +
  theme_bw(base_size = 13) +
  labs(
    title = "Relative abundance of top taxa by sample type",
    x = "Sample Type",
    y = "Relative abundance",
    fill = "Genus"
  )


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

PS_rare <- rarefy_even_depth(
  PS_16S_no_control,
  rngseed = 123,
  sample.size = 29600,
  replace = FALSE,
  verbose = FALSE
)

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


### BLAST INVESTIGATION 

tax_df <- as.data.frame(tax_table(PS_16S))
  
chloroplast_asvs <- rownames(tax_df)[
  apply(tax_df, 1, function(x)
    any(grepl("chloroplast", x, ignore.case = TRUE))
    )
  ]
  
length(chloroplast_asvs)

chloro_seqs <- rep_seqs[names(rep_seqs) %in% chloroplast_asvs]

length(chloro_seqs)

library(Biostrings)

writeXStringSet(
  chloro_seqs,
  filepath = "chloroplast_ASVs.fasta",
  format = "fasta"
)

  

#Alpha by time 

PS_sponge_rare

PS_sponge_rarefied



###############################################################################
# Alpha diversity (Observed ASVs) — Control vs Inoculated
#
# Purpose:
#   - Compare observed ASV richness between Control and Inoculated samples
#     across timepoints
#   - Uses rarefied data to control for sequencing depth
#   - Displays Wilcoxon rank-sum test results on the plot
#
# Input:
#   - PS_Object (phyloseq object of interest)
###############################################################################

library(phyloseq)
library(tidyverse)
library(ggpubr)

###############################################################################
# Alpha diversity (Observed ASVs) — Sponge samples
# Control vs Inoculated across timepoints
#
# Input:
#   - ps_sponge_rarefied_no_control (already rarefied)
###############################################################################

library(phyloseq)
library(tidyverse)
library(ggpubr)

sample_data(ps_sponge_rarefied_no_control)


# 1) Calculate Observed ASVs
alpha_df <- estimate_richness(
  ps_sponge_rarefied_no_control,
  measures = "Observed"
) %>%
  tibble::rownames_to_column("Sample")

# 2) Add metadata
meta_df <- data.frame(sample_data(ps_sponge_rarefied_no_control)) %>%
  tibble::rownames_to_column("Sample")

alpha_df <- left_join(alpha_df, meta_df, by = "Sample")

# 3) Keep only Control vs Inoculated and target timepoints
alpha_df <- alpha_df %>%
  filter(
    Tank_Type %in% c("Control", "Inoculated", "Existing_tank"),
    Timepoint %in% c("Pre_inoculation", "1_week_post_inoculation", "Harvest")
  )

# (Optional) set timepoint order for plotting
alpha_df$Timepoint <- factor(
  alpha_df$Timepoint,
  levels = c("Pre_inoculation", "1_week_post_inoculation", "Harvest")
)

TANK_COLORS <- c(
  "Existing_tank" = "#7A7A7A",
  "Control"       = "#F28E2B",
  "Inoculated"    = "#4E79A7"
)
ggplot(alpha_df, aes(x = Tank_Type, y = Observed, fill = Tank_Type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  facet_wrap(~ Timepoint, scales = "free_x") +
  
  # Control vs Inoculated (all timepoints)
  stat_compare_means(
    method = "wilcox.test",
    comparisons = list(c("Control", "Inoculated")),
    label = "p.format"
  ) +
  
  scale_fill_manual(values = TANK_COLORS) +
  
  # Existing tank comparisons (Harvest only, staggered)
  stat_compare_means(
    data = subset(alpha_df, Timepoint == "Harvest"),
    method = "wilcox.test",
    comparisons = list(
      c("Existing_tank", "Control"),
      c("Existing_tank", "Inoculated")
    ),
    label = "p.format",
    label.y = c(
      max(alpha_df$Observed) * 1.05,
      max(alpha_df$Observed) * 1.15
    )
  ) +
  
  theme_bw(base_size = 13) +
  labs(
    x = NULL,
    y = "Observed ASVs",
    title = "Observed ASV richness in sponge samples"
  ) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey90")
  )


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
  ps_sponge_rarefied,
  function(x) as.numeric(x > 0)
)

# 2) Compute Jaccard distance
jaccard_dist <- phyloseq::distance(
  ps_sponge_rarefied_no_control,
  method = "jaccard"
)

# 3) PCoA ordination
ord_jaccard <- ordinate(
  ps_sponge_rarefied_no_control,
  method = "PCoA",
  distance = jaccard_dist
)

# 4) Extract ordination scores + metadata
ord_df <- plot_ordination(
  PS_sponge_rarefied,
  ord_jaccard,
  justDF = TRUE
)

# 5) PCoA plot with confidence ellipses
p <- ggplot(
  ord_df,
  aes(x = Axis.1, y = Axis.2, color = Timepoint)
) +
  geom_point(size = 3, alpha = 0.9) +
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
    color = "Tank_Type"
  )

p

meta_df <- data.frame(sample_data(ps_sponge_rarefied_no_control))

adonis_jaccard <- adonis2(
  jaccard_dist ~ Timepoint,
  data = meta_df,
  permutations = 999
)

adonis_jaccard

anova(betadisper(jaccard_dist, meta_df$Timepoint))

###############################################################################
# Beta diversity (Jaccard presence/absence) — Tank_Type comparison
###############################################################################

library(phyloseq)
library(tidyverse)
library(vegan)

# Presence / absence transform
ps_pa <- transform_sample_counts(
  Ps_sponge_rarefied_no_control,
  function(x) as.numeric(x > 0)
)

# Jaccard distance
jaccard_dist <- phyloseq::distance(
  ps_sponge_rarefied_no_control,
  method = "jaccard"
)

# PCoA ordination
ord_jaccard <- ordinate(
  ps_sponge_rarefied_no_control,
  method = "PCoA",
  distance = jaccard_dist
)

# Extract scores + metadata
ord_df <- plot_ordination(
  ps_sponge_rarefied_no_control,
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


meta_df <- data.frame(sample_data(ps_sponge_rarefied_no_control))

adonis_tank <- adonis2(
  jaccard_dist ~ Tank_Type,
  data = meta_df,
  permutations = 999
)

adonis_tank

anova(betadisper(jaccard_dist, meta_df$Tank_Type))

