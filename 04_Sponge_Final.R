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


##ADJUSTED TO RELFECT INDIVIDUAL SAMPLES 

# 1) Transform to relative abundance
PS_rel <- transform_sample_counts(
  PS_16S,
  function(x) x / sum(x)
)

# 2) Identify top 50 taxa globally (you commented top 30 but used 50)
top_taxa <- names(sort(taxa_sums(PS_rel), decreasing = TRUE))[1:50]

# 3) Prune phyloseq object to top taxa only
PS_rel_top <- prune_taxa(top_taxa, PS_rel)

# 4) Melt to long format
rel_df <- psmelt(PS_rel_top)

# OPTIONAL: make sure Genus NA is labeled
rel_df$Genus <- ifelse(is.na(rel_df$Genus), "Unclassified", rel_df$Genus)

# 5) Aggregate within *individual samples* (NOT Sample_Type)
rel_df_sum <- rel_df %>%
  group_by(Sample, Sample_Type, Genus) %>%   # <-- key change
  summarise(Abundance = sum(Abundance), .groups = "drop")

# 6) Re-normalize so each *sample* sums to 1
rel_df_sum <- rel_df_sum %>%
  group_by(Sample) %>%                       # <-- key change
  mutate(Abundance = Abundance / sum(Abundance)) %>%
  ungroup()

# 7) Colors
n_taxa <- length(unique(rel_df_sum$Genus))
distinct_colors <- palette36.colors(n_taxa)
names(distinct_colors) <- unique(rel_df_sum$Genus)

# 8) Plot (individual samples on x-axis)
ggplot(rel_df_sum, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_bar(
    stat = "identity",
    color = "black",
    linewidth = 0.15
  ) +
  facet_grid(
    ~ Sample_Type,
    scales = "free_x",
    space = "free_x"
  ) +
  scale_fill_manual(values = distinct_colors) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom",
    legend.box = "horizontal"
  ) +
  labs(
    title = "Relative abundance of top taxa by sample",
    x = NULL,
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
  PS_sponge_rarefied,
  method = "jaccard"
)

# 3) PCoA ordination
ord_jaccard <- ordinate(
  PS_sponge_rarefied,
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
  aes(x = Axis.1, y = Axis.2, color = Tank_Type)
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


library(ggrepel)

p <- ggplot(
  ord_df,
  aes(x = Axis.1, y = Axis.2, color = Tank_Type)
) +
  geom_point(size = 3, alpha = 0.9) +
  geom_text_repel(
    aes(label = Sample),      # change if your column is named differently
    size = 3,
    max.overlaps = Inf,
    show.legend = FALSE
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

p

colnames(ord_df)

ord_df

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


###############################################################################
# Phylum-level relative abundance with Streptomyces highlighted
#
# Purpose:
#   - Visualize the relative abundance of Actinobacteriota across timepoints
#   - Highlight Streptomyces as a distinct slice
#   - Group remaining Actinobacteriota together
#   - Collapse all other taxa (including NA taxonomy) into "Other"
#   - Uses NON-rarefied data (appropriate for relative abundance)
#
# Input:
#   - PS_16S_no_control (phyloseq object)
###############################################################################

library(phyloseq)
library(tidyverse)

# 1) Transform counts to relative abundance (per sample)
ps_rel <- transform_sample_counts(
  PS_sponge_no_control,
  function(x) x / sum(x)
)

# 2) Melt to long format
rel_df <- psmelt(ps_rel)

# 3) Define plotting groups
rel_df <- rel_df %>%
  mutate(
    Plot_Group = case_when(
      Phylum == "Actinobacteriota" & Genus == "Streptomyces" ~ "Streptomyces",
      Phylum == "Actinobacteriota" ~ "Other_Actinobacteriota",
      TRUE ~ "Other"
    )
  )

# 4) Keep target timepoints and enforce order
rel_df <- rel_df %>%
  filter(Timepoint %in% c(
    "Pre_inoculation",
    "1_week_post_inoculation",
    "Harvest"
  )) %>%
  mutate(
    Timepoint = factor(
      Timepoint,
      levels = c(
        "Pre_inoculation",
        "1_week_post_inoculation",
        "Harvest"
      )
    )
  )

# 5) Aggregate and re-normalize within each timepoint
pie_df <- rel_df %>%
  group_by(Timepoint, Plot_Group) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  group_by(Timepoint) %>%
  mutate(Abundance = Abundance / sum(Abundance)) %>%
  ungroup()

# 6) Plot pie charts
ggplot(pie_df, aes(x = "", y = Abundance, fill = Plot_Group)) +
  geom_col(width = 1, color = "black") +
  coord_polar(theta = "y") +
  facet_wrap(~ Timepoint) +
  scale_fill_manual(
    values = c(
      "Streptomyces" = "#0072B2",
      "Other_Actinobacteriota" = "#D55E00",
      "Other" = "grey80"
    )
  ) +
  theme_void(base_size = 13) +
  labs(
    title = "Relative abundance of Actinobacteriota and Streptomyces across timepoints",
    fill = NULL
  )


###############################################################################
# Phylum-level relative abundance across timepoints (pie charts)
#
# Purpose:
#   - Visualize overall microbial community composition at the phylum level
#   - One pie chart per timepoint (Pre, Post, Harvest)
#   - Uses NON-rarefied data (appropriate for relative abundance)
#
# Input:
#   - PS_16S_no_control (phyloseq object)
###############################################################################

library(phyloseq)
library(tidyverse)
library(colorspace)

# 1) Transform counts to relative abundance (per sample)
ps_rel <- transform_sample_counts(
  Ps_sponge_no_euk_no_control,
  function(x) x / sum(x)
)

# 2) Melt to long format
rel_df <- psmelt(ps_rel)

# 3) Keep target timepoints and enforce order
rel_df <- rel_df %>%
  filter(Timepoint %in% c(
    "Pre_inoculation",
    "1_week_post_inoculation",
    "Harvest"
  )) %>%
  mutate(
    Timepoint = factor(
      Timepoint,
      levels = c(
        "Pre_inoculation",
        "1_week_post_inoculation",
        "Harvest"
      )
    ),
    Phylum = ifelse(is.na(Phylum), "Unclassified", Phylum)
  )

# 4) Aggregate and re-normalize within each timepoint
pie_df <- rel_df %>%
  group_by(Timepoint, Phylum) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  group_by(Timepoint) %>%
  mutate(Abundance = Abundance / sum(Abundance)) %>%
  ungroup()

# 5) Generate scalable qualitative colors
library(colorspace)

# Order phyla by overall abundance (most abundant first)
phylum_order <- pie_df %>%
  group_by(Phylum) %>%
  summarise(total_abundance = sum(Abundance), .groups = "drop") %>%
  arrange(desc(total_abundance)) %>%
  pull(Phylum)

pie_df$Phylum <- factor(pie_df$Phylum, levels = phylum_order)

# Generate maximally distinct colors
phylum_colors <- setNames(
  qualitative_hcl(
    n = length(phylum_order),
    palette = "Glasbey"
  ),
  phylum_order
)

# 6) Plot pie charts
ggplot(pie_df, aes(x = "", y = Abundance, fill = Phylum)) +
  geom_col(width = 1, color = "black", linewidth = 0.2) +
  coord_polar(theta = "y") +
  facet_wrap(~ Timepoint) +
  scale_fill_manual(values = phylum_colors) +
  theme_void(base_size = 13) +
  labs(
    title = "Phylum-level composition of sponge microbiomes across timepoints",
    fill = "Phylum"
  )


viewtax_table(PS_sponge_no_control)

sample_data(PS_sponge_no_control)

sample_sums(PS_sponge_no_control)

PS_sponge_no_control <- subset_samples(PS_sponge,Sample_Number != "73")

tax_df_PS_sponge_no_control <- as.data.frame(tax_table(PS_sponge_no_control))

Ps_sponge_no_euk_no_control <- subset_taxa(PS_sponge_no_control,Kingdom !="d__Eukaryota" )


library(phyloseq)
library(tidyverse)
library(colorspace)

# 1) Transform counts to relative abundance (per sample)
ps_rel <- transform_sample_counts(
  Ps_sponge_no_euk_no_control,
  function(x) x / sum(x)
)

# 2) Melt to long format
rel_df <- psmelt(ps_rel)

# 3) Keep target timepoints and clean taxonomy
rel_df <- rel_df %>%
  filter(Timepoint %in% c(
    "Pre_inoculation",
    "1_week_post_inoculation",
    "Harvest"
  )) %>%
  mutate(
    Timepoint = factor(
      Timepoint,
      levels = c(
        "Pre_inoculation",
        "1_week_post_inoculation",
        "Harvest"
      )
    ),
    Phylum = ifelse(is.na(Phylum), "Unclassified", Phylum)
  )

# 4) Aggregate by Timepoint × Tank_Type × Phylum
pie_df <- rel_df %>%
  group_by(Timepoint, Tank_Type, Phylum) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  group_by(Timepoint, Tank_Type) %>%
  mutate(Abundance = Abundance / sum(Abundance)) %>%
  ungroup()

# 5) Order phyla globally for consistent colors
phylum_order <- pie_df %>%
  group_by(Phylum) %>%
  summarise(total_abundance = sum(Abundance), .groups = "drop") %>%
  arrange(desc(total_abundance)) %>%
  pull(Phylum)

pie_df <- pie_df %>%
  mutate(
    Tank_Type = factor(
      Tank_Type,
      levels = c(
        "Control",
        "Inoculated",
        "Existing_tank"
      )
    )
  )

phylum_colors <- setNames(
  qualitative_hcl(
    n = length(phylum_order),
    palette = "Dynamic"   # use Dynamic; Glasbey not valid here
  ),
  phylum_order
)

library(Polychrome)

phylum_colors <- setNames(
  createPalette(
    length(levels(pie_df$Phylum)),
    c("#000000", "#FFFFFF")
  ),
  levels(pie_df$Phylum)
)


ggplot(pie_df, aes(x = "", y = Abundance, fill = Phylum)) +
  geom_col(width = 1.1, color = "black", linewidth = 0.2) +
  coord_polar(theta = "y") +
  facet_grid(
    Tank_Type ~ Timepoint,
    drop = TRUE,
    switch = "y"
  ) +
  scale_fill_manual(values = phylum_colors) +
  theme_void(base_size = 13) +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    
    strip.text.x = element_text(face = "bold", size = 13),
    strip.text.y = element_text(face = "bold", size = 13),
    strip.placement = "outside",
    
    panel.spacing.x = unit(2.2, "lines"),
    panel.spacing.y = unit(1.8, "lines"),
    
    plot.margin = margin(20, 20, 20, 20)
  ) +
  labs(fill = "Phylum")


ggsave(
  "phylum_pies_by_time_and_tank.png",
  width = 14,
  height = 8,
  dpi = 600
)


###############################################################################
# Genus-level relative abundance (Top 20 genera + Other)
#
# Purpose:
#   - Visualize genus-level community composition using NON-rarefied data
#   - Explicitly label unclassified genera (NA → "Unclassified")
#   - Keep only the top 20 most abundant genera
#   - Collapse all remaining genera into "Other"
#   - Maintain consistent genus colors across plots
#
# Input:
#   - PS_sponge_no_euk (phyloseq object, Eukaryota removed)
###############################################################################

library(phyloseq)
library(tidyverse)
library(colorspace)

# 1) Transform to relative abundance (per sample)
ps_rel <- transform_sample_counts(
  Ps_sponge_no_euk_no_control,
  function(x) x / sum(x)
)

# 2) Melt to long format
rel_df <- psmelt(ps_rel)

# 3) Explicitly label unclassified genus
rel_df <- rel_df %>%
  mutate(
    Genus = ifelse(is.na(Genus), "Unclassified", Genus)
  )

# 4) Keep only Control vs Inoculated
rel_df <- rel_df %>%
  filter(Tank_Type %in% c("Control", "Inoculated"))

# 5) Identify top 20 genera globally
top20_genera <- rel_df %>%
  group_by(Genus) %>%
  summarise(total_abundance = sum(Abundance), .groups = "drop") %>%
  arrange(desc(total_abundance)) %>%
  slice_head(n = 20) %>%
  pull(Genus)

# 6) Collapse non-top genera into "Other"
rel_df <- rel_df %>%
  mutate(
    Genus_plot = ifelse(Genus %in% top20_genera, Genus, "Other")
  )

# 7) Aggregate by Timepoint, Tank_Type, Genus
rel_df_sum <- rel_df %>%
  group_by(Timepoint, Tank_Type, Genus_plot) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

# 8) Re-normalize within each Timepoint × Tank_Type
rel_df_sum <- rel_df_sum %>%
  group_by(Timepoint, Tank_Type) %>%
  mutate(Abundance = Abundance / sum(Abundance)) %>%
  ungroup()

# 9) Lock genus order and colors globally
genus_levels <- rel_df_sum %>%
  group_by(Genus_plot) %>%
  summarise(total_abundance = sum(Abundance), .groups = "drop") %>%
  arrange(desc(total_abundance)) %>%
  pull(Genus_plot)

rel_df_sum$Genus_plot <- factor(rel_df_sum$Genus_plot, levels = genus_levels)

genus_colors <- c(
  setNames(
    qualitative_hcl(
      n = length(genus_levels) - 1,
      palette = "Glasbey"
    ),
    setdiff(genus_levels, "Other")
  ),
  "Other" = "grey80"
)

# ---- rel_df_sum is now ready for plotting ----



ggplot(
  filter(rel_df_sum, Timepoint == "Pre_inoculation"),
  aes(x = Tank_Type, y = Abundance, fill = Genus_plot)
) +
  geom_col(color = "black", linewidth = 0.15) +
  scale_fill_manual(values = genus_colors) +
  theme_bw(base_size = 13) +
  labs(
    title = "Pre-inoculation",
    x = NULL,
    y = "Relative abundance",
    fill = "Genus"
  )


ggplot(
  filter(rel_df_sum, Timepoint == "1_week_post_inoculation"),
  aes(x = Tank_Type, y = Abundance, fill = Genus_plot)
) +
  geom_col(color = "black", linewidth = 0.15) +
  scale_fill_manual(values = genus_colors) +
  theme_bw(base_size = 13) +
  labs(
    title = "Pre-inoculation",
    x = NULL,
    y = "Relative abundance",
    fill = "Genus"
  )


ggplot(
  filter(rel_df_sum, Timepoint == "Harvest"),
  aes(x = Tank_Type, y = Abundance, fill = Genus_plot)
) +
  geom_col(color = "black", linewidth = 0.15) +
  scale_fill_manual(values = genus_colors) +
  theme_bw(base_size = 13) +
  labs(
    title = "Pre-inoculation",
    x = NULL,
    y = "Relative abundance",
    fill = "Genus"
  )

length(unique(rel_df$Genus))

rel_df %>%
  group_by(Genus) %>%
  summarise(total = sum(Abundance)) %>%
  arrange(desc(total)) %>%
  head(30)



###############################################################################
# Genus-level relative abundance (Top 15 genera + Other)
#
# Purpose:
#   - Visualize genus-level community composition using NON-rarefied data
#   - Keep top 15 genera globally
#   - Collapse remaining genera into "Other"
#   - Maintain consistent genus colors across plots
#   - Generate three plots:
#       1) Pre-inoculation (Control vs Inoculated)
#       2) Post-inoculation (Control vs Inoculated)
#       3) Harvest (Control vs Inoculated vs Existing_tank)
#
# Input:
#   - PS_sponge_no_euk (phyloseq object, Eukaryota removed)
###############################################################################

library(phyloseq)
library(tidyverse)
library(colorspace)
library(Polychrome)


# 1) Transform to relative abundance (per sample)
ps_rel <- transform_sample_counts(
  Ps_sponge_no_euk_no_control,
  function(x) x / sum(x)
)

# 2) Melt to long format
rel_df <- psmelt(ps_rel)

#3) Label unclassified genera explicitly
rel_df <- rel_df %>%
  mutate(Genus = ifelse(is.na(Genus), "Unclassified", Genus))

# 4) Identify top 15 genera globally
top15_genera <- rel_df %>%
  group_by(Genus) %>%
  summarise(total_abundance = sum(Abundance), .groups = "drop") %>%
  arrange(desc(total_abundance)) %>%
  slice_head(n = 30) %>%
  pull(Genus)

# 5) Collapse non-top genera into "Other"
rel_df <- rel_df %>%
  mutate(
    Genus_plot = ifelse(Genus %in% top15_genera, Genus, "Other")
  )

# 6) Aggregate by Timepoint, Tank_Type, Genus
rel_df_sum <- rel_df %>%
  group_by(Timepoint, Tank_Type, Genus_plot) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

# 7) Re-normalize within each Timepoint × Tank_Type
rel_df_sum <- rel_df_sum %>%
  group_by(Timepoint, Tank_Type) %>%
  mutate(Abundance = Abundance / sum(Abundance)) %>%
  ungroup()

# 8) Lock genus order and colors globally
genus_levels <- rel_df_sum %>%
  group_by(Genus_plot) %>%
  summarise(total_abundance = sum(Abundance), .groups = "drop") %>%
  arrange(desc(total_abundance)) %>%
  pull(Genus_plot)

rel_df_sum$Genus_plot <- factor(rel_df_sum$Genus_plot, levels = genus_levels)

# Create Glasbey-style colors for genera
genus_colors <- c(
  setNames(
    createPalette(
      length(genus_levels) - 1,
      c("#000000", "#FFFFFF")  # seed colors to maximize contrast
    ),
    setdiff(genus_levels, "Other")
  ),
  "Other" = "grey80"
)


# 9) Plotting function
plot_timepoint <- function(tp, title_text, tanks) {
  ggplot(
    filter(rel_df_sum, Timepoint == tp, Tank_Type %in% tanks),
    aes(x = Tank_Type, y = Abundance, fill = Genus_plot)
  ) +
    geom_col(color = "black", linewidth = 0.15) +
    scale_fill_manual(values = genus_colors) +
    theme_bw(base_size = 13) +
    labs(
      title = title_text,
      x = NULL,
      y = "Relative abundance",
      fill = "Genus"
    )
}

# 10) Generate plots
p_pre <- plot_timepoint(
  "Pre_inoculation",
  "Pre-inoculation",
  c("Control", "Inoculated")
)

p_post <- plot_timepoint(
  "1_week_post_inoculation",
  "Post-inoculation",
  c("Control", "Inoculated")
)

p_harvest <- plot_timepoint(
  "Harvest",
  "Harvest",
  c("Control", "Inoculated", "Existing_tank")
)

# Display plots
p_pre
p_post
p_harvest


library(patchwork)

(p_pre + p_post + p_harvest) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

sort(unique(rel_df_sum$Genus_plot))
sort(names(genus_colors))

setdiff(
  unique(rel_df_sum$Genus_plot),
  names(genus_colors)
)


library(phyloseq)
library(tidyverse)
library(colorspace)
library(Polychrome)

# 1) Relative abundance per sample
ps_rel <- transform_sample_counts(
  Ps_sponge_no_euk_no_control,
  function(x) x / sum(x)
)

# 2) Melt
rel_df <- psmelt(ps_rel) %>%
  mutate(Genus = ifelse(is.na(Genus), "Unclassified", Genus))

# 3) Top 30 genera globally
top_genera <- rel_df %>%
  group_by(Genus) %>%
  summarise(total_abundance = sum(Abundance), .groups = "drop") %>%
  arrange(desc(total_abundance)) %>%
  slice_head(n = 30) %>%
  pull(Genus)

# 4) Collapse others
rel_df <- rel_df %>%
  mutate(Genus_plot = ifelse(Genus %in% top_genera, Genus, "Other"))

# 5) Keep INDIVIDUAL samples
rel_df_sum <- rel_df %>%
  group_by(Sample, Timepoint, Tank_Type, Genus_plot) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  group_by(Sample) %>%
  mutate(Abundance = Abundance / sum(Abundance)) %>%
  ungroup()

# 6) Order Tank_Type
rel_df_sum <- rel_df_sum %>%
  mutate(
    Tank_Type = factor(Tank_Type, levels = c("Control", "Inoculated", "Existing_tank"))
  )

# 7) Order Timepoint (THIS IS THE FIX)
rel_df_sum <- rel_df_sum %>%
  mutate(
    Timepoint = factor(
      Timepoint,
      levels = c("Pre_inoculation", "1_week_post_inoculation", "Harvest")
    )
  )

# 8) Order samples within each timepoint by Tank_Type
sample_levels <- rel_df_sum %>%
  distinct(Sample, Timepoint, Tank_Type) %>%
  arrange(Timepoint, Tank_Type, Sample) %>%
  pull(Sample)

rel_df_sum <- rel_df_sum %>%
  mutate(Sample = factor(Sample, levels = sample_levels))

# 9) Genus order + colors
genus_levels <- rel_df_sum %>%
  group_by(Genus_plot) %>%
  summarise(total_abundance = sum(Abundance), .groups = "drop") %>%
  arrange(desc(total_abundance)) %>%
  pull(Genus_plot)

rel_df_sum$Genus_plot <- factor(rel_df_sum$Genus_plot, levels = genus_levels)

genus_colors <- c(
  setNames(
    createPalette(length(genus_levels) - 1, c("#000000", "#FFFFFF")),
    setdiff(genus_levels, "Other")
  ),
  "Other" = "grey80"
)

# 10) Plot: individual samples, Tank_Type labels, facetted by Timepoint
ggplot(
  rel_df_sum,
  aes(x = Sample, y = Abundance, fill = Genus_plot)
) +
  geom_col(color = "black", linewidth = 0.15) +
  facet_wrap(
    ~ Timepoint,
    nrow = 1,
    scales = "free_x"
  ) +
  scale_x_discrete(
    labels = function(x) {
      rel_df_sum$Tank_Type[match(x, rel_df_sum$Sample)]
    }
  ) +
  scale_fill_manual(values = genus_colors) +
  theme_bw(base_size = 13) +
  theme(
    # ---- axes ----
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 9
    ),
    
    # ---- legend (VERTICAL) ----
    legend.position = "bottom",
    legend.box = "vertical",
    legend.key.size = unit(0.45, "cm"),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 12),
    legend.spacing.y = unit(0.3, "cm"),
    
    # ---- titles ----
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5)
  ) +
  guides(
    fill = guide_legend(
      ncol = 1,        # <<< THIS makes labels vertical
      byrow = TRUE
    )
  ) +
  labs(
    title = "Timepoint",
    subtitle = "Genus-level relative abundance by individual sample",
    x = "Tank Type",
    y = "Relative abundance",
    fill = "Genus"
  )


ggplot(
  rel_df_sum,
  aes(x = Sample, y = Abundance, fill = Genus_plot)
) +
  geom_col(color = "black", linewidth = 0.15) +
  facet_wrap(
    ~ Timepoint,
    nrow = 1,
    scales = "free_x"
  ) +
  scale_x_discrete(
    labels = function(x) {
      rel_df_sum$Tank_Type[match(x, rel_df_sum$Sample)]
    }
  ) +
  scale_fill_manual(values = genus_colors) +
  theme_bw(base_size = 13) +
  theme(
    # ---- axes ----
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 9
    ),
    
    # ---- legend (VERTICAL) ----
    legend.position = "bottom",
    legend.box = "vertical",
    legend.key.size = unit(0.45, "cm"),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 12),
    legend.spacing.y = unit(0.3, "cm"),
    
    # ---- titles ----
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5)
  ) +
  guides(
    fill = guide_legend(
      ncol = 1,        # <<< THIS makes labels vertical
      byrow = TRUE
    )
  ) +
  labs(
    title = "Timepoint",
    subtitle = "Genus-level relative abundance by individual sample",
    x = "Tank Type",
    y = "Relative abundance",
    fill = "Genus"
  )




ggplot(
  rel_df_sum,
  aes(x = Sample, y = Abundance, fill = Genus_plot)
) +
  geom_col(color = "black", linewidth = 0.15) +
  facet_wrap(
    ~ Timepoint,
    nrow = 1,
    scales = "free_x"
  ) +
  scale_x_discrete(
    labels = function(x) {
      rel_df_sum$Tank_Type[match(x, rel_df_sum$Sample)]
    }
  ) +
  scale_fill_manual(values = genus_colors) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 9
    ),
    
    # ---- LEGEND: make genus labels larger ----
    legend.position = "bottom",
    legend.box = "vertical",
    legend.key.size = unit(0.55, "cm"),   # slightly larger color boxes
    legend.text = element_text(size = 13),# <<< BIGGER GENUS LABELS
    legend.title = element_text(size = 14),
    legend.spacing.y = unit(0.35, "cm"),
    
    # ---- titles ----
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    plot.subtitle = element_text(size = 13, hjust = 0.5)
  ) +
  guides(
    fill = guide_legend(
      ncol = 1,        # keep legend vertical
      byrow = TRUE
    )
  ) +
  labs(
    title = "Timepoint",
    subtitle = "Genus-level relative abundance by individual sample",
    x = "Tank Type",
    y = "Relative abundance",
    fill = "Genus"
  )

ggplot(
  rel_df_sum,
  aes(x = Sample, y = Abundance, fill = Genus_plot)
) +
  geom_col(color = "black", linewidth = 0.15) +
  facet_wrap(
    ~ Timepoint,
    nrow = 1,
    scales = "free_x"
  ) +
  scale_x_discrete(
    labels = function(x) {
      rel_df_sum$Tank_Type[match(x, rel_df_sum$Sample)]
    }
  ) +
  scale_fill_manual(values = genus_colors) +
  theme_bw(base_size = 13) +
  theme(
    # X-axis labels vertical
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 9
    ),
    
    # FORCE legend to be horizontal
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.key.size = unit(0.55, "cm"),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 14),
    legend.spacing.x = unit(0.6, "cm"),
    legend.spacing.y = unit(0.2, "cm"),
    
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    plot.subtitle = element_text(size = 13, hjust = 0.5)
  ) +
  guides(
    fill = guide_legend(
      nrow = 1,      # <<< THIS IS THE KEY LINE
      byrow = TRUE
    )
  ) +
  labs(
    title = "Timepoint",
    subtitle = "Genus-level relative abundance by individual sample",
    x = "Tank Type",
    y = "Relative abundance",
    fill = "Genus"
  )

ggplot(
  rel_df_sum,
  aes(x = Sample, y = Abundance, fill = Genus_plot)
) +
  geom_col(color = "black", linewidth = 0.15) +
  facet_wrap(
    ~ Timepoint,
    nrow = 1,
    scales = "free_x"
  ) +
  scale_x_discrete(
    labels = function(x) {
      rel_df_sum$Tank_Type[match(x, rel_df_sum$Sample)]
    }
  ) +
  scale_fill_manual(values = genus_colors) +
  theme_bw(base_size = 13) +
  theme(
    # X-axis labels vertical
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 9
    ),
    
    # FORCE legend to be horizontal
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.key.size = unit(0.55, "cm"),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 14),
    legend.spacing.x = unit(0.6, "cm"),
    legend.spacing.y = unit(0.2, "cm"),
    
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    plot.subtitle = element_text(size = 13, hjust = 0.5)
  ) +
  guides(
    fill = guide_legend(
      nrow = 1,      # <<< THIS IS THE KEY LINE
      byrow = TRUE
    )
  ) +
  labs(
    title = "Timepoint",
    subtitle = "Genus-level relative abundance by individual sample",
    x = "Tank Type",
    y = "Relative abundance",
    fill = "Genus"
  )


ggplot(
  rel_df_sum,
  aes(x = Sample, y = Abundance, fill = Genus_plot)
) +
  geom_col(color = "black", linewidth = 0.15) +
  facet_wrap(
    ~ Timepoint,
    nrow = 1,
    scales = "free_x"
  ) +
  scale_x_discrete(
    labels = function(x) {
      rel_df_sum$Tank_Type[match(x, rel_df_sum$Sample)]
    }
  ) +
  scale_fill_manual(values = genus_colors) +
  theme_bw(base_size = 13) +
  theme(
    # --- make plot visually smaller ---
    panel.spacing = unit(0.8, "lines"),
    plot.margin = margin(t = 8, r = 8, b = 2, l = 8),
    
    # --- x-axis labels ---
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 9
    ),
    
    # --- legend formatting ---
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.key.size = unit(0.5, "cm"),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 14),
    legend.spacing.x = unit(0.4, "cm"),
    legend.spacing.y = unit(0.3, "cm"),
    
    # --- titles ---
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    plot.subtitle = element_text(size = 13, hjust = 0.5)
  ) +
  guides(
    fill = guide_legend(
      nrow = 3,        # <<< STACK legend into multiple rows
      byrow = TRUE
    )
  ) +
  labs(
    title = "Timepoint",
    subtitle = "Genus-level relative abundance by individual sample",
    x = "Tank Type",
    y = "Relative abundance",
    fill = "Genus"
  )

ggsave(
  filename = "genus_relative_abundance.png",
  plot = last_plot(),
  width = 20,
  height = 10,
  units = "in",
  dpi = 300
)

ggsave("figure.pdf", width = 25, height = 10)




###############################################################################
# Genus-level relative abundance using cumulative abundance cutoff (90%)
#
# Purpose:
#   - Visualize genus-level community composition using NON-rarefied data
#   - Retain as many genera as needed to explain 90% of total abundance
#   - Collapse remaining low-abundance genera into "Other"
#   - Maintain consistent genus colors across plots
#
# Input:
#   - PS_sponge_no_euk (phyloseq object, Eukaryota removed)
###############################################################################

library(phyloseq)
library(tidyverse)
library(Polychrome)

# 1) Transform to relative abundance (per sample)
ps_rel <- transform_sample_counts(
  Ps_sponge_no_euk_no_control,
  function(x) x / sum(x)
)

# 2) Melt to long format
rel_df <- psmelt(ps_rel)

# 3) Label unclassified genera explicitly
rel_df <- rel_df %>%
  mutate(Genus = ifelse(is.na(Genus), "Unclassified", Genus))

# 4) Compute global genus abundance
genus_abundance <- rel_df %>%
  group_by(Genus) %>%
  summarise(total_abundance = sum(Abundance), .groups = "drop") %>%
  arrange(desc(total_abundance)) %>%
  mutate(
    cumulative_abundance = cumsum(total_abundance) / sum(total_abundance)
  )

# 5) Select genera up to 90% cumulative abundance
keep_genera <- genus_abundance %>%
  filter(cumulative_abundance <= 0.90) %>%
  pull(Genus)

# Ensure Unclassified is kept if abundant
keep_genera <- union(keep_genera, "Unclassified")

# 6) Collapse remaining genera into "Other"
rel_df <- rel_df %>%
  mutate(
    Genus_plot = ifelse(Genus %in% keep_genera, Genus, "Other")
  )

# 7) Aggregate by Timepoint and Tank_Type
rel_df_sum <- rel_df %>%
  group_by(Timepoint, Tank_Type, Genus_plot) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

# 8) Re-normalize within each Timepoint × Tank_Type
rel_df_sum <- rel_df_sum %>%
  group_by(Timepoint, Tank_Type) %>%
  mutate(Abundance = Abundance / sum(Abundance)) %>%
  ungroup()

# 9) Lock genus order and colors globally
genus_levels <- rel_df_sum %>%
  group_by(Genus_plot) %>%
  summarise(total_abundance = sum(Abundance), .groups = "drop") %>%
  arrange(desc(total_abundance)) %>%
  pull(Genus_plot)

rel_df_sum$Genus_plot <- factor(rel_df_sum$Genus_plot, levels = genus_levels)

genus_colors <- c(
  setNames(
    createPalette(
      length(genus_levels) - 1,
      c("#000000", "#FFFFFF")
    ),
    setdiff(genus_levels, "Other")
  ),
  "Other" = "grey80"
)

length(keep_genera)

keep_genera

###############################################################################
# Genus-level relative abundance with biologically meaningful "Other" groups
#
# Purpose:
#   - Visualize genus-level community composition using NON-rarefied data
#   - Avoid a single dominant "Other" category
#   - Split low-abundance genera into phylum-aware groups
#   - Maintain consistent colors across all plots
#
# Input:
#   - PS_sponge_no_euk (phyloseq object, Eukaryota removed)
###############################################################################

library(phyloseq)
library(tidyverse)
library(Polychrome)

# 1) Transform to relative abundance (per sample)
ps_rel <- transform_sample_counts(
  Ps_sponge_no_euk_no_control,
  function(x) x / sum(x)
)

# 2) Melt to long format
rel_df <- psmelt(ps_rel)

# 3) Label unclassified genus
rel_df <- rel_df %>%
  mutate(
    Genus = ifelse(is.na(Genus), "Unclassified", Genus)
  )

# 4) Define plotting groups (genus + phylum-aware Others)
rel_df <- rel_df %>%
  mutate(
    Genus_plot = case_when(
      Genus == "Unclassified" ~ "Unclassified",
      Phylum == "Actinobacteriota" ~ paste0("Other Actinobacteriota"),
      Phylum == "Proteobacteria" ~ paste0("Other Proteobacteria"),
      Phylum == "Firmicutes" ~ paste0("Other Firmicutes"),
      TRUE ~ "Other (remaining phyla)"
    )
  )

# Keep named genera separate
named_genera <- unique(rel_df$Genus)
rel_df$Genus_plot[rel_df$Genus %in% named_genera] <- rel_df$Genus[rel_df$Genus %in% named_genera]

# 5) Aggregate
rel_df_sum <- rel_df %>%
  group_by(Timepoint, Tank_Type, Genus_plot) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

# 6) Re-normalize within Timepoint × Tank_Type
rel_df_sum <- rel_df_sum %>%
  group_by(Timepoint, Tank_Type) %>%
  mutate(Abundance = Abundance / sum(Abundance)) %>%
  ungroup()

# 7) Lock order and colors
genus_levels <- rel_df_sum %>%
  group_by(Genus_plot) %>%
  summarise(total_abundance = sum(Abundance), .groups = "drop") %>%
  arrange(desc(total_abundance)) %>%
  pull(Genus_plot)

rel_df_sum$Genus_plot <- factor(rel_df_sum$Genus_plot, levels = genus_levels)

genus_colors <- setNames(
  createPalette(length(genus_levels), c("#000000", "#FFFFFF")),
  genus_levels
)


###############################################################################
# Relative abundance of selected genera across timepoints
#
# Purpose:
#   - Highlight specific genera of interest
#   - Compare Control vs Inoculated across Pre and Post timepoints
#   - Include Existing tank at Harvest
#   - Collapse all non-target genera into "Other"
#   - Ensure relative abundances sum to 100%
#
# Input:
#   - PS_sponge_no_euk (phyloseq object, NON-rarefied)
###############################################################################

library(phyloseq)
library(tidyverse)
library(Polychrome)

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
  Ps_sponge_no_euk_no_control,
  function(x) x / sum(x)
)

# 2) Melt to long format
rel_df <- psmelt(ps_rel)

# 3) Handle unclassified genus
rel_df <- rel_df %>%
  mutate(Genus = ifelse(is.na(Genus), "Unclassified", Genus))

# 4) Define plotting groups: genera of interest vs Other
rel_df <- rel_df %>%
  mutate(
    Genus_plot = ifelse(
      Genus %in% genera_of_interest,
      Genus,
      "Other"
    )
  )

# 5) Aggregate by Timepoint, Tank_Type, Genus_plot
rel_df_sum <- rel_df %>%
  group_by(Timepoint, Tank_Type, Genus_plot) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

# 6) Re-normalize within each Timepoint × Tank_Type
rel_df_sum <- rel_df_sum %>%
  group_by(Timepoint, Tank_Type) %>%
  mutate(Abundance = Abundance / sum(Abundance)) %>%
  ungroup()

# 7) Lock factor order and colors
genus_levels <- c(genera_of_interest, "Other")

rel_df_sum$Genus_plot <- factor(rel_df_sum$Genus_plot, levels = genus_levels)

genus_colors <- c(
  setNames(
    createPalette(length(genera_of_interest), c("#000000", "#FFFFFF")),
    genera_of_interest
  ),
  "Other" = "grey80"
)

# 8) Plotting function
plot_timepoint <- function(tp, title_text, tanks) {
  ggplot(
    filter(rel_df_sum, Timepoint == tp, Tank_Type %in% tanks),
    aes(x = Tank_Type, y = Abundance, fill = Genus_plot)
  ) +
    geom_col(color = "black", linewidth = 0.15) +
    scale_fill_manual(values = genus_colors) +
    theme_bw(base_size = 13) +
    labs(
      title = title_text,
      x = NULL,
      y = "Relative abundance",
      fill = "Genus"
    )
}

# 9) Generate plots
p_pre <- plot_timepoint(
  "Pre_inoculation",
  "Pre-inoculation",
  c("Control", "Inoculated")
)

p_post <- plot_timepoint(
  "1_week_post_inoculation",
  "Post-inoculation",
  c("Control", "Inoculated")
)

p_harvest <- plot_timepoint(
  "Harvest",
  "Harvest",
  c("Control", "Inoculated", "Existing_tank")
)

# Display plots
p_pre
p_post
p_harvest



###############################################################################
# Pie charts of selected genera across timepoints
#
# Purpose:
#   - Highlight selected genera of interest
#   - Compare Control vs Inoculated (Pre, Post)
#   - Include Existing tank at Harvest
#   - Collapse all non-target genera into "Other"
#   - Ensure relative abundances sum to 100%
#
# Input:
#   - PS_sponge_no_euk (phyloseq object, NON-rarefied)
###############################################################################

library(phyloseq)
library(tidyverse)
library(Polychrome)

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
  Ps_sponge_no_euk_no_control,
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

