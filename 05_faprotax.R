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
  "C:\\Users\\mhewapat\\Documents\\avery_Farms\\16S_full_run\\R_project\\16S_Avery_Full_Run\\faprotax_export\\taxonomy.tsv",
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE
)

library(phyloseq)

ps <- phyloseq(
  otu_table(as.matrix(otu_raw), taxa_are_rows = TRUE),
  tax_table(as.matrix(tax))
)

ps

sample_names(ps)
taxa_names(ps)[1:10]

head(colnames(otu_raw))

head(tax$Taxon)


head(meta)

colnames(meta)

sample_names(ps)

sample_data(ps) <- sample_data(meta)

ps


identical(rownames(sample_data(ps)), sample_names(ps))
colnames(sample_data(ps))



library(stringr)
library(phyloseq)
library(dplyr)


tax_df <- as.data.frame(tax_table(ps))

tax_df$Genus <- str_extract(tax_df$Taxon, "g__[^;]+")
tax_df$Genus <- str_remove(tax_df$Genus, "g__")
tax_df$Genus[tax_df$Genus == ""] <- NA

tax_table(ps) <- tax_table(as.matrix(tax_df))

rank_names(ps)

head(tax_table(ps)[, c("Taxon", "Genus")])

colnames(tax_table(ps))

head(tax_table(ps)[, "Genus"])
table(is.na(tax_table(ps)[, "Genus"]))


ps_g <- subset_taxa(ps, !is.na(Genus) & Genus != "Unassigned")

ps_g

ps_genus <- tax_glom(ps_g, taxrank = "Genus")
ps_genus <- prune_taxa(taxa_sums(ps_genus) > 0, ps_genus)

ps_genus


faprotax_input <- as.data.frame(otu_table(ps_genus))
faprotax_input$Genus <- as.vector(tax_table(ps_genus)[, "Genus"])

# Put Genus as the first column
faprotax_input <- faprotax_input[, c("Genus", setdiff(colnames(faprotax_input), "Genus"))]

head(faprotax_input)
dim(faprotax_input)

colnames(faprotax_input)


faprotax_input <- faprotax_input[, -1]

write.table(
  faprotax_input,
  "faprotax_input.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)



head(readLines("faprotax_input.tsv", n = 5))

head(faprotax_input$Genus, 20)

length(unique(faprotax_input$Genus))


sum(faprotax_input$Genus %in% c(
  "Streptomyces",
  "Pseudomonas",
  "Bacillus",
  "Rhizobium",
  "Nitrosomonas",
  "Nitrospira"
))

rowSums(faprotax_input[faprotax_input$Genus %in% c("Pseudomonas", "Bacillus"), -1])


