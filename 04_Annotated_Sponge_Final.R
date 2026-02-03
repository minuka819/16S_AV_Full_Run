

### RELATIVE ABUNDANCE PLOTS 

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



library(dplyr)
library(ggplot2)

## 1) Transform to relative abundance
PS_rel <- transform_sample_counts(
  PS_16S,
  function(x) x / sum(x)
)

## 2) Identify top 50 taxa globally
top_taxa <- names(sort(taxa_sums(PS_rel), decreasing = TRUE))[1:50]

## 3) Prune phyloseq object to top taxa only
PS_rel_top <- prune_taxa(top_taxa, PS_rel)

## 4) Melt to long format
rel_df <- psmelt(PS_rel_top)

## Label unclassified genera
rel_df$Genus <- ifelse(is.na(rel_df$Genus), "Unclassified", rel_df$Genus)

## 5) Aggregate within *individual samples*
rel_df_sum <- rel_df %>%
  group_by(Sample, Sample_Type, Tank_Type, Timepoint, Genus) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

## 6) Re-normalize so each *sample* sums to 1
rel_df_sum <- rel_df_sum %>%
  group_by(Sample) %>%
  mutate(Abundance = Abundance / sum(Abundance)) %>%
  ungroup()

## 7) Create combined Tank + Timepoint label
rel_df_sum <- rel_df_sum %>%
  mutate(
    Sample_tag = paste(Tank_Type, Timepoint, sep = " | ")
  )


library(dplyr)

# Get x positions for tank annotations
tank_labels <- rel_df_sum %>%
  distinct(Sample, Sample_Type, Tank_Type) %>%
  filter(Sample_Type == "Sponge") %>%
  mutate(x_pos = as.numeric(factor(Sample))) %>%
  group_by(Tank_Type) %>%
  summarise(
    x = mean(x_pos),
    .groups = "drop"
  )


## 8) Colors
n_taxa <- length(unique(rel_df_sum$Genus))
distinct_colors <- palette36.colors(n_taxa)
names(distinct_colors) <- unique(rel_df_sum$Genus)

ggplot(rel_df_sum,
       aes(
         x = interaction(Tank_Type, Sample, drop = TRUE),
         y = Abundance,
         fill = Genus
       )) +
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
  scale_x_discrete(
    labels = function(x) sub("\\..*", "", x)  # keep only Tank_Type
  ) +
  scale_fill_manual(values = distinct_colors) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1,
      size = 9
    ),
    legend.position = "bottom",
    legend.box = "horizontal"
  ) +
  labs(
    title = "Relative abundance of top taxa by individual sample",
    x = "Tank Type",
    y = "Relative abundance",
    fill = "Genus"
  )

library(dplyr)
library(ggplot2)

library(dplyr)

rel_df_sum <- rel_df_sum %>%
  mutate(
    Tank_Type = case_when(
      Tank_Type %in% c("Existing", "Existing_tank", "Existing tank") ~ "Existing",
      Tank_Type %in% c("Control", "control") ~ "Control",
      Tank_Type %in% c("Inoculated", "inoculated") ~ "Inoculated",
      TRUE ~ Tank_Type
    )
  )



# 1) Define order
tank_order <- c("Control", "Inoculated", "Existing")

# 2) Ensure Tank_Type is NEVER NA and is ordered
rel_df_sum <- rel_df_sum %>%
  mutate(
    Tank_Type = ifelse(is.na(Tank_Type) | Tank_Type == "", "Unknown", as.character(Tank_Type)),
    Tank_Type = factor(Tank_Type, levels = c(tank_order, "Unknown"))
  )

# 3) Create a per-sample ordering table (unique rows per sample)
sample_order <- rel_df_sum %>%
  distinct(Sample, Sample_Type, Tank_Type) %>%
  arrange(Sample_Type, Tank_Type, Sample) %>%
  mutate(Sample = factor(Sample, levels = unique(Sample)))

# 4) Join ordering back (keeps all genera rows, but fixes bar order + labeling)
rel_df_sum <- rel_df_sum %>%
  left_join(sample_order, by = c("Sample", "Sample_Type", "Tank_Type"))

# 5) Plot: each bar is still an individual Sample, but the LABEL shown is Tank_Type
ggplot(rel_df_sum, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.15) +
  facet_grid(~ Sample_Type, scales = "free_x", space = "free_x") +
  scale_x_discrete(labels = function(x) {
    # label each sample by its Tank_Type (no collapsing)
    lbl <- sample_order$Tank_Type[match(x, as.character(sample_order$Sample))]
    as.character(lbl)
  }) +
  scale_fill_manual(values = distinct_colors) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 9),
    legend.position = "bottom",
    legend.box = "horizontal"
  ) +
  labs(
    title = "Relative abundance of top taxa by individual sample",
    x = "Tank Type (each bar is an individual sample)",
    y = "Relative abundance",
    fill = "Genus"
  )

library(dplyr)

tank_order <- c("Control", "Inoculated", "Existing_tank")

sample_order <- rel_df_sum %>%
  distinct(Sample, Sample_Type, Tank_Type) %>%
  mutate(
    Tank_Type = factor(Tank_Type, levels = tank_order)
  ) %>%
  arrange(Sample_Type, Tank_Type, Sample) %>%   # <-- ordering logic
  mutate(
    Sample_ordered = factor(Sample, levels = Sample)
  ) %>%
  select(Sample, Sample_ordered)

rel_df_sum <- rel_df_sum %>%
  left_join(sample_order, by = "Sample")

ggplot(rel_df_sum,
       aes(x = Sample_ordered, y = Abundance, fill = Genus)) +
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
  scale_x_discrete(
    labels = function(x) {
      rel_df_sum$Tank_Type[match(x, rel_df_sum$Sample_ordered)]
    }
  ) +
  scale_fill_manual(values = distinct_colors) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1,
      size = 9
    ),
    legend.position = "bottom",
    legend.box = "horizontal"
  ) +
  labs(
    title = "Relative abundance of top taxa by individual sample",
    x = "Tank Type (bars ordered: Control → Inoculated → Existing)",
    y = "Relative abundance",
    fill = "Genus"
  )


library(dplyr)
library(ggplot2)

# --------------------------------------------
# 1) CLEAN + FORCE Tank_Type (THIS IS THE FIX)
# --------------------------------------------
rel_df_sum <- rel_df_sum %>%
  mutate(
    Tank_Type = as.character(Tank_Type),
    Tank_Type = trimws(Tank_Type)
  )

# Explicitly map the three known types
rel_df_sum <- rel_df_sum %>%
  mutate(
    Tank_Type = case_when(
      Tank_Type == "Control"    ~ "Control",
      Tank_Type == "Inoculated" ~ "Inoculated",
      is.na(Tank_Type)          ~ "Existing_tank",
      TRUE                      ~ Tank_Type
    )
  )

# Lock the order (AFTER mapping)
rel_df_sum <- rel_df_sum %>%
  mutate(
    Tank_Type = factor(
      Tank_Type,
      levels = c("Control", "Inoculated", "Existing_tank")
    )
  )

# --------------------------------------------
# 2) CREATE ORDERED SAMPLE FACTOR
# --------------------------------------------
sample_levels <- rel_df_sum %>%
  distinct(Sample, Sample_Type, Tank_Type) %>%
  arrange(Sample_Type, Tank_Type, Sample) %>%
  pull(Sample)

rel_df_sum <- rel_df_sum %>%
  mutate(
    Sample_ordered = factor(Sample, levels = sample_levels)
  )

# --------------------------------------------
# 3) PLOT (each bar = sample, labeled by Tank_Type)
# --------------------------------------------
ggplot(
  rel_df_sum,
  aes(x = Sample_ordered, y = Abundance, fill = Genus)
) +
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
  scale_x_discrete(
    labels = function(x) {
      rel_df_sum$Tank_Type[match(x, rel_df_sum$Sample)]
    }
  ) +
  scale_fill_manual(values = distinct_colors) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1,
      size = 9
    ),
    legend.position = "bottom",
    legend.box = "horizontal"
  ) +
  labs(
    title = "Relative abundance of top taxa by individual sample",
    x = "Tank Type (Control → Inoculated → Existing tank)",
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
    color = "Timepoint"
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



#bray curtis



# 1) Compute Bray–Curtis distance (no rarefaction needed)
bray_dist <- phyloseq::distance(
  PS_16S,          # or PS_abund_filt if you choose to filter
  method = "bray"
)

# 2) PCoA ordination
ord_bray <- ordinate(
  PS_16S,
  method = "PCoA",
  distance = bray_dist
)

# 3) Extract ordination scores + metadata
ord_df <- plot_ordination(
  PS_16S,
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


library(phyloseq)
library(ggplot2)

bray_dist <- phyloseq::distance(PS_16S, method = "bray")

ord_bray <- ordinate(
  PS_16S,
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
