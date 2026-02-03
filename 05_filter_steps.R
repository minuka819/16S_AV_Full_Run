

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

PS_abund_filt <- prune_taxa(asv_percent_total >= 0.01, PS_16S)

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

# Read depth per sample (after abundance filtering)
sample_reads_filt <- sample_sums(PS_abund_filt)

# Summary stats
summary(sample_reads_filt)

# Quantiles (same ones you used before)
quantile(sample_reads_filt, probs = c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1))

par(mfrow = c(1, 2))

hist(sample_sums(PS_16S),
     breaks = 30,
     main = "Before filtering",
     xlab = "Reads per sample")

hist(sample_reads_filt,
     breaks = 30,
     main = "After 0.01% filtering",
     xlab = "Reads per sample")

par(mfrow = c(1, 1))


library(dplyr)
library(ggplot2)

# Create read depth table
read_depth_df <- data.frame(
  SampleID = names(sample_sums(PS_16S)),
  Reads_Before = sample_sums(PS_16S),
  Reads_After  = sample_sums(PS_abund_filt)
)

# Percent retained per sample
read_depth_df$Percent_Retained <- 
  (read_depth_df$Reads_After / read_depth_df$Reads_Before) * 100

# Inspect
head(read_depth_df)
summary(read_depth_df$Percent_Retained)


ggplot(read_depth_df, aes(y = reorder(SampleID, Reads_Before))) +
  geom_segment(aes(x = Reads_Before, xend = Reads_After, yend = SampleID),
               color = "grey60") +
  geom_point(aes(x = Reads_Before), color = "black", size = 2) +
  geom_point(aes(x = Reads_After), color = "steelblue", size = 2) +
  theme_bw() +
  labs(
    title = "Per-sample read depth before and after ASV abundance filtering",
    x = "Read count",
    y = "Sample"
  )


library(dplyr)
library(tidyr)
library(ggplot2)

read_depth_df <- data.frame(
  SampleID = names(sample_sums(PS_16S)),
  Before = sample_sums(PS_16S),
  After  = sample_sums(PS_abund_filt)
)

# Convert to long format for ggplot
read_depth_long <- read_depth_df %>%
  pivot_longer(
    cols = c(Before, After),
    names_to = "Stage",
    values_to = "Reads"
  )


ggplot(read_depth_long, aes(x = SampleID, y = Reads, fill = Stage)) +
  geom_col(position = "dodge") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  labs(
    title = "Per-sample read depth before and after ASV abundance filtering",
    x = "Sample",
    y = "Read depth",
    fill = "Filtering stage"
  )





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
  PS_abund_filt,
  method = "PCoA",
  distance = jaccard_dist
)

# 4) Extract ordination scores + metadata
ord_df <- plot_ordination(
  PS_abund_filt,
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


