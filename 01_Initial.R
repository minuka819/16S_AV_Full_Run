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


PS_16S

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

# Sort by Reads (ascending)
read_depth_all <- read_depth_all[order(read_depth_all$Reads), ]

# Force ggplot to respect this order
read_depth_all$Sample <- factor(
  read_depth_all$Sample,
  levels = read_depth_all$Sample
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
final_plot <- plot_read_depth_before + arrow_plot + plot_read_depth_after +
  plot_layout(widths = c(1, 0.12, 1))

final_plot


