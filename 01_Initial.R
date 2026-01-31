### --- CHAPTER 1 --- ###

### --- INITIAL QC AND INVESTIGATIONS --- ###

# Helper function: install and load required packages automatically
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

rep_seqs <- read_qza(
  "C:\\Users\\mhewapat\\Documents\\avery_Farms\\16S_full_run\\results\\artifacts\\repseqs_16S.qza"
)$data

refseq(PS_16S) <- rep_seqs

# 1️⃣ Extract OTU table
otu <- as(otu_table(PS_16S), "matrix")
if (taxa_are_rows(PS_16S)) otu <- t(otu)

#Get metadata 
meta <- as(sample_data(PS_16S), "data.frame") %>%
  tibble::rownames_to_column("SampleID")

get_curve <- function(counts, sample_id) {
  depths <- seq(100, sum(counts), length.out = 50)   # avoid 0 to prevent rarefy issues
  observed <- sapply(depths, function(d) {
    rarefy(rbind(counts), sample = d, se = FALSE)
  })
  data.frame(SampleID = sample_id, Depth = depths, Observed = observed)
}

rare_df <- map_dfr(rownames(otu), ~get_curve(otu[.x, ], .x))

rare_df <- left_join(rare_df, meta, by = "SampleID")

ggplot(rare_df, aes(x = Depth, y = Observed, group = SampleID, color = Sample_Type)) +
  geom_line(alpha = 0.8, linewidth = 1) +
  facet_wrap(~ Sample_Type, scales = "free_y") +
  theme_bw(base_size = 13) +
  labs(
    title = "Rarefaction Curves by Sample Type",
    x = "Sequencing depth (reads)",
    y = "Observed Features",
    color = "Sample Type"
  ) +
 scale_color_manual(values = c("Leaf" = "#028A0F", "Sponge" = "#02A3D3")) +
 scale_x_continuous(labels = scales::comma)


PS_16S_Sponge <- subset_samples(PS_16S,Sample_Type == "Sponge")
PS_16S_Leaf <- subset_samples(PS_16S,Sample_Type == "Leaf")

meta_leaf <- sample_data(PS_16S_Leaf)

meta_sponge <- sample_data(PS_16S_Sponge)

PS_16S_Existing <- subset_samples(PS_16S,Tank_Type == "Existing_tank")

view(sample_data(PS_16S_Existing), "data.frame")%>%
  tibble::rownames_to_column("SampleID")

PS_16S_Control <- subset_samples(PS_16S,Tank_Type == "Control")

view(sample_data(PS_16S_Control), "data.frame")%>%
  tibble::rownames_to_column("SampleID")

PS_16S_Inoculated <- subset_samples(PS_16S,Tank_Type == "Inoculated")

view(sample_data(PS_16S_Control), "data.frame")%>%
  tibble::rownames_to_column("SampleID")


read_depths <- data.frame(
  Sample       = sample_names(PS_16S_Sponge),
  Tank_Type    = sample_data(PS_16S_Sponge)$Tank_Type,
  Timepoint    = sample_data(PS_16S_Sponge)$Timepoint,
  Reads        = sample_sums(PS_16S_Sponge)
)

read_depths <- read_depths %>%
  arrange(Reads) %>%
  mutate(Sample = factor(Sample, levels = Sample))

ggplot(read_depths, aes(x = Sample, y = Reads, fill = Tank_Type)) +
  geom_col() +
  facet_wrap(~ Timepoint, scales = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(y = "Read Depth", x = "Sample", fill = "Tank Type")


#### ------------------ [03] QC Filtering ------------------ ####
# 1) Keep only bacterial ASVs (adjust "Kingdom" if your column name differs)

min_depth <- min(sample_sums(ps_mocks_fungi))   # should be ~35k+

PS_16S_no_control <- subset_samples(PS_16S, Sample_Number != "73")

view(sample_data(PS_16S_no_control), "data.frame")%>%
  tibble::rownames_to_column("SampleID")

PS_16S_no_control <- prune_taxa(taxa_sums(PS_16S_no_control) > 0, PS_16S_no_control)


depths <- sample_sums(PS_16S_no_control)
min_depth <- min(depths)
min_depth

library(phyloseq)
library(ggplot2)

# Build read depth dataframe
read_depth_all <- data.frame(
  Sample = sample_names(PS_16S),
  Reads = sample_sums(PS_16S),
  Sample_Type = sample_data(PS_16S)$Sample_Type,
  Timepoint = sample_data(PS_16S)$Timepoint
)


# Plot

ggplot(read_depth_all, aes(x = Sample, y = Reads, fill = Sample_Type)) +
  geom_col() +
  scale_fill_manual(
    values = c(
      "Leaf" = "#028A0F",
      "Sponge" = "#02A3D3"
    )
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(
    title = "Read depth per sample",
    x = "Sample",
    y = "Read depth",
    fill = "Sample Type"
  )



# Subset to sponge + leaf
ps_qc <- subset_samples(PS_16S, Sample_Type %in% c("Sponge", "Leaf"))
ps_qc <- prune_taxa(taxa_sums(ps_qc) > 0, ps_qc)

# Build dataframe
read_depth_df <- data.frame(
  Reads = sample_sums(ps_qc),
  Sample_Type = sample_data(ps_qc)$Sample_Type,
  Timepoint = sample_data(ps_qc)$Timepoint
)

# Ensure timepoint order (adjust if needed)
read_depth_df$Timepoint <- factor(
  read_depth_df$Timepoint,
  levels = c("Pre_inoculation", "1_week_post_inoculation", "Harvest")
)

# Plot
ggplot(read_depth_df, aes(x = Timepoint, y = Reads, fill = Sample_Type)) +
  geom_col(position = "dodge") +
  scale_fill_manual(
    values = c(
      "Leaf" = "#028A0F",
      "Sponge" = "#02A3D3"
    )
  ) +
  theme_bw() +
  labs(
    title = "Read depth across timepoints (pre-filtering)",
    x = "Timepoint",
    y = "Read depth",
    fill = "Sample Type"
  )+
  scale_y_continuous(labels = scales::comma)

library(dplyr)

read_depth_summary <- read_depth_df %>%
  group_by(Timepoint, Sample_Type) %>%
  summarise(
    mean_reads = mean(Reads),
    sd_reads   = sd(Reads),
    .groups = "drop"
  )


ggplot(read_depth_summary, aes(x = Timepoint, y = mean_reads, fill = Sample_Type)) +
  geom_col(position = position_dodge(width = 0.7)) +
  geom_errorbar(
    aes(
      ymin = mean_reads - sd_reads,
      ymax = mean_reads + sd_reads
    ),
    width = 0.2,
    position = position_dodge(width = 0.7)
  ) +
  scale_fill_manual(
    values = c(
      "Leaf" = "#028A0F",
      "Sponge" = "#02A3D3"
    )
  ) +
  theme_bw(base_size = 13) +
  labs(
    title = "Read depth across timepoints (pre-filtering)",
    x = "Timepoint",
    y = "Read depth",
    fill = "Sample Type"
  ) +
  scale_y_continuous(labels = scales::comma)



ggplot(read_depth_df, aes(x = Timepoint, y = Reads, color = Sample_Type)) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.4,
    position = position_dodge(width = 0.7)
  ) +
  geom_jitter(
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.7
    ),
    alpha = 0.8,
    size = 2
  ) +
  scale_color_manual(
    values = c(
      "Leaf" = "#028A0F",
      "Sponge" = "#02A3D3"
    )
  ) +
  theme_bw(base_size = 13) +
  labs(
    title = "Per-sample read depth by timepoint",
    x = "Timepoint",
    y = "Read depth",
    color = "Sample Type"
  ) +
  scale_y_continuous(labels = scales::comma)


library(phyloseq)
library(ggplot2)

PS_16S_no_control

# Build read depth dataframe
read_depth_PS_no_control <- data.frame(
  Sample = sample_names(PS_16S_no_control),
  Reads = sample_sums(PS_16S_no_control),
  Sample_Type = sample_data(PS_16S_no_control)$Sample_Type,
  Timepoint = sample_data(PS_16S_no_control)$Timepoint
)

# Plot

ggplot(read_depth_PS_no_control, aes(x = Sample, y = Reads, fill = Sample_Type)) +
  geom_col() +
  scale_fill_manual(
    values = c(
      "Leaf" = "#028A0F",
      "Sponge" = "#02A3D3"
    )
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(
    title = "Read depth per sample",
    x = "Sample",
    y = "Read depth",
    fill = "Sample Type"
  )

depth_df <- data.frame(
  Sample = sample_names(PS_16S_no_control),
  Reads = sample_sums(PS_16S_no_control),
  Sample_Type = sample_data(PS_16S_no_control)$Sample_Type
)

tapply(depths, sample_data(PS_16S_no_control)$Sample_Type, summary)
tapply(depths, sample_data(PS_16S_no_control)$Sample_Type, min)

29603


set.seed(123)
RARE_DEPTH <- 29600  # example

PS_rare <- rarefy_even_depth(
  PS_16S_no_control,
  sample.size = RARE_DEPTH,
  rngseed = 123,
  replace = FALSE,
  verbose = FALSE
)

PS_rare

library(phyloseq)
library(ggplot2)
library(patchwork)

# ---- Helper: make a read depth df ----
make_depth_df <- function(ps) {
  df <- data.frame(
    Sample = sample_names(ps),
    Reads = sample_sums(ps),
    Sample_Type = sample_data(ps)$Sample_Type
  )
  df <- df[order(df$Reads), ]
  df$Sample <- factor(df$Sample, levels = df$Sample)
  df
}

# ---- BEFORE (pre-rarefaction) ----
ymax <- max(sample_sums(PS_16S_no_control))
ymax

depth_before <- make_depth_df(PS_16S_no_control)

plot_read_depth_before <- ggplot(depth_before, aes(x = Sample, y = Reads, fill = Sample_Type)) +
  geom_col() +
  scale_fill_manual(values = c("Leaf" = "#028A0F", "Sponge" = "#02A3D3")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Before rarefaction", x = "Sample", y = "Read depth", fill = "Sample Type")

# ---- AFTER (post-rarefaction) ----
depth_after <- make_depth_df(PS_rare)

plot_read_depth_after <- ggplot(depth_after, aes(x = Sample, y = Reads, fill = Sample_Type)) +
  geom_col() +
  scale_fill_manual(values = c("Leaf" = "#028A0F", "Sponge" = "#02A3D3")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "After rarefaction (29,600 reads/sample)", x = "Sample", y = "Read depth", fill = "Sample Type")

# ---- Arrow plot (your style) ----
arrow_plot <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = "→", size = 20) +
  theme_void()

# ---- Final side-by-side ----

plot_read_depth_before <- plot_read_depth_before +
  coord_cartesian(ylim = c(0, ymax))

plot_read_depth_after <- plot_read_depth_after +
  coord_cartesian(ylim = c(0, ymax))

final_plot <- plot_read_depth_before +scale_y_continuous(labels = scales::comma) + arrow_plot + plot_read_depth_after +
  plot_layout(widths = c(1, 0.12, 1)) + scale_y_continuous(labels = scales::comma)


final_plot


###-Alpha Diversity-###



library(phyloseq)
library(ggplot2)

alpha_df <- estimate_richness(PS_rare, measures = c("Observed", "Shannon"))
alpha_df$Sample_Type <- sample_data(PS_rare)$Sample_Type

p_obs <- wilcox.test(Observed ~ Sample_Type, data = alpha_df)$p.value
p_sha <- wilcox.test(Shannon  ~ Sample_Type, data = alpha_df)$p.value

# Observed ASVs
ggplot(alpha_df, aes(x = Sample_Type, y = Observed, fill = Sample_Type)) +
  geom_boxplot() +
  theme_bw() +
  labs(title = "Alpha diversity (Observed ASVs) - rarefied", x = NULL, y = "Observed ASVs") +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("Leaf" = "#028A0F", "Sponge" = "#02A3D3"))+

# Shannon
ggplot(alpha_df, aes(x = Sample_Type, y = Shannon, fill = Sample_Type)) +
  geom_boxplot() +
  theme_bw() +
  labs(title = "Alpha diversity (Shannon) - rarefied", x = NULL, y = "Shannon") +
  theme(legend.position = "none")+
  scale_fill_manual(values = c("Leaf" = "#028A0F", "Sponge" = "#02A3D3"))

COLORS <- c(
  "Leaf"   = "#028A0F",
  "Sponge" = "#02A3D3"
)


p1 <- ggplot(alpha_df, aes(x = Sample_Type, y = Observed, fill = Sample_Type)) +
  geom_boxplot() +
  scale_fill_manual(values = COLORS) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(title = "Observed ASVs", x = NULL, y = "Observed") +
  annotate(
    "text",
    x = 1.5,
    y = max(alpha_df$Observed) * 1.05,
    label = paste0("Wilcoxon p = ", formatC(p_obs, format = "e", digits = 2))
  )

p2 <- ggplot(alpha_df, aes(x = Sample_Type, y = Shannon, fill = Sample_Type)) +
  geom_boxplot() +
  scale_fill_manual(values = COLORS) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(title = "Shannon diversity", x = NULL, y = "Shannon") +
  annotate(
    "text",
    x = 1.5,
    y = max(alpha_df$Shannon) * 1.05,
    label = paste0("Wilcoxon p = ", formatC(p_sha, format = "e", digits = 2))
  )

p1 + p2

library(phyloseq)
library(ggplot2)

COLORS <- c("Leaf" = "#028A0F", "Sponge" = "#02A3D3")

# Use non-rarefied for composition plots
ps <- PS_16S_no_control

# Keep only Leaf + Sponge
ps <- subset_samples(ps, Sample_Type %in% c("Leaf", "Sponge"))
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

# (Recommended) remove host organelles if present in taxonomy
# If this errors due to rank names, run: rank_names(ps)
ps <- subset_taxa(ps, Order != "Chloroplast" & Family != "Mitochondria")
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

# Agglomerate to Genus (change to "Phylum" if you prefer)
ps_g <- tax_glom(ps, taxrank = "Genus", NArm = FALSE)

# Convert to relative abundance
ps_rel <- transform_sample_counts(ps_g, function(x) x / sum(x))

# Identify Top 30 taxa across all samples
top30 <- names(sort(taxa_sums(ps_rel), decreasing = TRUE))[1:30]

# Prune to Top 30 and collapse the rest into "Other"
ps_top <- prune_taxa(top30, ps_rel)

# Melt to long format for ggplot
df <- psmelt(ps_top)

# If Genus is NA, label as "Unassigned"
df$Genus <- as.character(df$Genus)
df$Genus[is.na(df$Genus) | df$Genus == ""] <- "Unassigned"

# Plot: bars per sample, faceted by Sample_Type (clean comparison)
ggplot(df, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_col() +
  facet_wrap(~ Sample_Type, scales = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(
    title = "Top 30 genera (relative abundance)",
    x = "Sample",
    y = "Relative abundance",
    fill = "Genus"
  )


library(dplyr)

df_sum <- df %>%
  group_by(Sample_Type, Genus) %>%
  summarise(MeanAbundance = mean(Abundance), .groups = "drop")

ggplot(df_sum, aes(x = Sample_Type, y = MeanAbundance, fill = Genus)) +
  geom_col() +
  scale_x_discrete(limits = c("Leaf", "Sponge")) +
  theme_bw() +
  labs(
    title = "Top 30 genera (mean relative abundance by sample type)",
    x = NULL,
    y = "Mean relative abundance",
    fill = "Genus"
  )

# ps_bact = after removing chloroplast + mito, before rarefaction
ps_bact <- subset_samples(PS_16S_no_control, Sample_Type %in% c("Leaf","Sponge"))
ps_bact <- prune_taxa(taxa_sums(ps_bact) > 0, ps_bact)
ps_bact <- subset_taxa(ps_bact, Order != "Chloroplast" & Family != "Mitochondria")
ps_bact <- prune_taxa(taxa_sums(ps_bact) > 0, ps_bact)

bact_depth <- data.frame(
  Sample = sample_names(ps_bact),
  Bacterial_Reads = sample_sums(ps_bact),
  Sample_Type = sample_data(ps_bact)$Sample_Type
)

bact_depth

ggplot(bact_depth, aes(x = Sample_Type, y = Bacterial_Reads, fill = Sample_Type)) +
  stat_summary(fun = mean, geom = "bar", width = 0.6) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  scale_fill_manual(values = COLORS) +
  theme_bw() +
  labs(
    title = "Total bacterial read depth after host removal",
    x = NULL,
    y = "Mean bacterial reads"
  ) +
  theme(legend.position = "none")

### - PIE CHARTS - ###

library(phyloseq)
library(dplyr)
library(tidyr)
library(ggplot2)

# ---------- Colors ----------
PIE_COLORS <- c("Host" = "grey70", "Bacterial" = "steelblue")

# ---------- Identify host taxa (chloroplast + mitochondria) ----------
is_chloro <- tax_table(PS_16S_no_control)[, "Order"] == "Chloroplast"
is_mito   <- tax_table(PS_16S_no_control)[, "Family"] == "Mitochondria"
host_ids  <- taxa_names(PS_16S_no_control)[(is_chloro | is_mito)]

# ---------- Per-sample host vs bacterial reads ----------
host_bact_df <- data.frame(
  Sample      = sample_names(PS_16S_no_control),
  Sample_Type = sample_data(PS_16S_no_control)$Sample_Type,
  Host_Reads  = sample_sums(prune_taxa(host_ids, PS_16S_no_control)),
  Total_Reads = sample_sums(PS_16S_no_control)
) %>%
  mutate(
    Bacterial_Reads = Total_Reads - Host_Reads
  )

# ---------- Mean fractions + mean reads by Sample_Type (for pie charts) ----------
pie_df <- host_bact_df %>%
  group_by(Sample_Type) %>%
  summarise(
    Host = mean(Host_Reads, na.rm = TRUE),
    Bacterial = mean(Bacterial_Reads, na.rm = TRUE),
    Total = mean(Total_Reads, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c("Host", "Bacterial"),
               names_to = "Category",
               values_to = "Reads") %>%
  mutate(
    Fraction = Reads / Total,
    Label = paste0(
      Category, "\n",
      format(round(Reads), big.mark = ","), " reads\n(",
      round(Fraction * 100, 1), "%)"
    )
  )

# ---------- Plot ----------
ggplot(pie_df, aes(x = "", y = Fraction, fill = Category)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  facet_wrap(~ Sample_Type) +
  geom_text(aes(label = Label), position = position_stack(vjust = 0.5), size = 4) +
  scale_fill_manual(values = PIE_COLORS) +
  theme_void() +
  labs(
    title = "Host vs bacterial reads by sample type",
    fill = NULL
  )



