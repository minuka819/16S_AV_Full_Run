# Load libraries
library(phyloseq)
library(stringr)
library(dplyr)

# ---- 1. Read feature table (FIXES V2–V74 ISSUE) ----
otu_raw <- read.table(
  "G:/Projects/AGC/metagenomics/16S_Full_Run/faprotax_export/feature-table.tsv",
  sep = "\t",
  header = TRUE,
  skip = 1,
  comment.char = "",
  check.names = FALSE
)

# Fix the first column name
colnames(otu_raw)[1] <- "ASV_ID"

# Move ASV_ID to rownames
rownames(otu_raw) <- otu_raw$ASV_ID
otu_raw$ASV_ID <- NULL

colnames(otu_raw)

# ---- 2. Read taxonomy table ----
tax <- read.table(
  "G:\\Projects\\AGC\\metagenomics\\16S_Full_Run\\faprotax_export\\taxonomy.tsv",
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE
)

# ---- 3. Build phyloseq object ----
ps <- phyloseq(
  otu_table(as.matrix(otu), taxa_are_rows = TRUE),
  tax_table(as.matrix(tax))
)

ps

# ---- 4. Extract Genus from Taxon ----
tax_df <- as.data.frame(tax_table(ps))

tax_df$Genus <- str_extract(tax_df$Taxon, "g__[^;]+")
tax_df$Genus <- str_remove(tax_df$Genus, "g__")
tax_df$Genus[tax_df$Genus == ""] <- NA

tax_table(ps) <- tax_table(as.matrix(tax_df))

# ---- 5. Remove taxa without usable genus ----
ps <- subset_taxa(
  ps,
  !is.na(Genus) &
    Genus != "Unassigned" &
    Genus != "Bacteria"
)

# ---- 6. Collapse ASVs → Genus ----
ps_genus <- tax_glom(ps, taxrank = "Genus")
ps_genus <- prune_taxa(taxa_sums(ps_genus) > 0, ps_genus)

# ---- 7. Create FAPROTAX input data frame ----
faprotax_input <- as.data.frame(otu_table(ps_genus))
faprotax_input$Genus <- as.vector(tax_table(ps_genus)[, "Genus"])

# Put Genus first
faprotax_input <- faprotax_input %>%
  select(Genus, everything())

faprotax_input

# ---- 8. View result ----
View(faprotax_input)
head(faprotax_input)
dim(faprotax_input)
