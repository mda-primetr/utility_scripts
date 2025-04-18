#!/usr/bin/env Rscript

# Process_16S_biom.R
#
# Usage:
#   Rscript Process_16S_biom.R /path/to/input.biom Project_name
#
# Reads a BIOM file and constructs a phyloseq object. Saves the phyloseq object
# as an RDS file for later use.
# Use Rbiom package rbiom_1.0.3
# Rest of the packages should be ok

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript Process_16S_biom.R <input_biom> <project_name>\n")
  quit(status = 1)
}

biom_file <- args[1]
project_name <- args[2]

# Load required packages
suppressPackageStartupMessages({
  library(rbiom)
  library(vegan)
  library(phyloseq)
  library(Biostrings)
  library(tidyr)
  library(dplyr)
  library(magrittr)
  library(tibble)
})



# By default this will use the sample with the lowest reads to rarefy all the samples
RNGkind(sample.kind = "Rounding")
set.seed(711)


# Read_biom file
d1 <- rbiom::read.biom(biom_file)


# function to create a phyloseq object from a BIOM file ----
create_phyloseq_16S <- function(d1_biom) {
  # Read BIOM file. This uses the biomformat package.


  # If your BIOM object does not store sequences in d1$sequences,
  # you'll need to adapt how you retrieve them:
  if (!is.null(d1_biom$sequences)) {
    seqs <- Biostrings::DNAStringSet(d1_biom$sequences)
  } else {
    seqs <- NULL
  }

  # Extract OTU (counts) matrix
  d1_OTU_mat <- as.matrix(rbiom::counts(d1_biom))

  # Extract taxonomy
  # Often, the first column can be an OTU or feature ID. Adjust indexing if needed.
  d1_TAX_mat <- as.matrix(taxonomy(d1_biom))[, -1]

  # Build phyloseq objects
  OTU <- otu_table(d1_OTU_mat, taxa_are_rows = TRUE)
  TAX <- tax_table(d1_TAX_mat)

  physeq <- phyloseq(OTU, TAX, seqs)

  # If a phylogenetic tree is available, add it
  if (!is.null(d1_biom$phylogeny)) {
    phy_tree(physeq) <- d1_biom$phylogeny
  }



  return(physeq)
}

#  Create the phyloseq object
physeq <- create_phyloseq_16S(d1)

# Save the phyloseq object as an RDS file
saveRDS(physeq, file = paste0(project_name, ".", format(Sys.Date(), "%Y.%m.%d"), ".", "phyloseq.rds"))



# Output Psmelt object
df_melted <- psmelt(physeq)
df_seqs <- data.frame(physeq@refseq) %>%
  tibble::rownames_to_column("ASVID") %>%
  rename("sequence" = 2) %>%
  mutate(sequence = as.character(sequence)) %>%
  rowwise() %>%
  mutate(sequence_hash = digest::digest(sequence, algo = "xxh3_128", serialize = FALSE)) %>%
  ungroup()


# wide format
df_melted %>%
  tidyr::pivot_wider(names_from = "Sample", values_from = "Abundance") %>%
  dplyr::rename(ASVID = OTU) %>%
  inner_join(df_seqs, by = "ASVID") %>%
  write.table(., file = paste0(project_name, ".", format(Sys.Date(), "%Y.%m.%d"), ".", "phyloseq_melted_wide.txt"), col.names = NA, sep = "\t")


df_melted %>%
  dplyr::rename(ASVID = OTU) %>%
  write.table(., file = paste0(project_name, ".", format(Sys.Date(), "%Y.%m.%d"), ".", ".phyloseq_melted_long.txt"), col.names = NA, sep = "\t")


# Melted tables with relative abundance
psmelt(physeq %>% transform_sample_counts(.,function(x) (x / sum(x))*100)) %>%
  dplyr::rename(ASVID = OTU) %>%
  write.table(., file = paste0(project_name, ".", format(Sys.Date(), "%Y.%m.%d"), ".", "relative_phyloseq_melted_long.txt"), col.names = NA, sep = "\t")


psmelt(physeq %>% transform_sample_counts(.,function(x) (x / sum(x))*100)) %>%
  dplyr::rename(ASVID = OTU) %>%
  tidyr::pivot_wider(names_from = "Sample", values_from = "Abundance") %>%
  write.table(., file = paste0(project_name, ".", format(Sys.Date(), "%Y.%m.%d"), ".", "relative_phyloseq_melted_wide.txt"), col.names = NA, sep = "\t")






# Prep data for alpha diversity
d1_OTU_mat <- as.matrix(rbiom::counts(d1))

# function  Alpha diversity ----
output_alpha_div <- function(OTU_mat, project = project_name, is_rarefied = FALSE) {
  # Determine the suffix based on whether it's rarefied or not
  suffix <- if (is_rarefied) "rarefied." else ""


  write.table(
    data.frame(
      cbind(
        number_of_classified_reads = colSums(OTU_mat),
        Observed = colSums(OTU_mat > 0),
        InvSimpson = vegan::diversity(t(OTU_mat), index = "invsimpson"),
        Shannon = vegan::diversity(t(OTU_mat), index = "shannon")
      )
    ),
    file = paste0(project_name, ".", format(Sys.Date(), "%Y.%m.%d"), ".", suffix, "no_of_classified_reads_and_alpha_diversity.txt"),
    col.names = NA, sep = "\t"
  )
}


# function Beta diversity ----

output_beta_div <- function(OTU_mat, project = project_name, is_rarefied = FALSE) {
  # Determine the suffix based on whether it's rarefied or not
  suffix <- if (is_rarefied) "rarefied." else ""

  # Bray distance as is
  write.table(as.matrix(vegan::vegdist(t(OTU_mat), method = "bray")),
    file = paste0(project, ".", format(Sys.Date(), "%Y.%m.%d"), ".", suffix, "bray_counts_dist.txt"), col.names = NA, sep = "\t"
  )

  # Bray distance on relative abundances
  write.table(as.matrix(vegan::vegdist(apply(t(OTU_mat), 2, function(x) (x / sum(x))), method = "bray")),
    file = paste0(project, ".", format(Sys.Date(), "%Y.%m.%d"), ".", suffix, "bray_relative_dist.txt"), col.names = NA, sep = "\t"
  )

  # Jaccard distance
  write.table(as.matrix(vegan::vegdist(t(OTU_mat), method = "jaccard")),
    file = paste0(project, ".", format(Sys.Date(), "%Y.%m.%d"), ".", suffix, "jaccard_dist.txt"), col.names = NA, sep = "\t"
  )

}




# Without Rarefaction ----

# Output Alpha diversity ----
output_alpha_div(d1_OTU_mat, project = project_name, is_rarefied = FALSE)


# Output beta diversity ----
output_beta_div(d1_OTU_mat, project = project_name, is_rarefied = FALSE)


# Beta diversity using unifrac distance ----
# Since UniFrac distance has not been implemented in the vegan package, we have to use Phyloseq package to derive it 
# Also note that it requires rooted tree and in this case phyloseq will randomly select one of the ASV as root

write.table(as.matrix(phyloseq::distance(physeq, method = "wunifrac")),
   file = paste0(project_name, ".", format(Sys.Date(), "%Y.%m.%d"), ".", "weighted_unifrac_dist.txt"), col.names = NA, sep = "\t"
 )





# Apply rarefaction ----

final_phyloseq_rarefied_16s <- rarefy_even_depth(physeq, rngseed = FALSE, verbose = FALSE)
# Save the phyloseq object as an RDS file
saveRDS(final_phyloseq_rarefied_16s, file = paste0(project_name, ".", format(Sys.Date(), "%Y.%m.%d"), ".", "phyloseq.rarefied.rds"))


# Output Psmelt object
df_melted_rare <- psmelt(final_phyloseq_rarefied_16s)

df_melted_rare %>%
  tidyr::pivot_wider(names_from = "Sample", values_from = "Abundance") %>%
  dplyr::rename(ASVID = OTU) %>%
  write.table(., file = paste0(project_name, ".", format(Sys.Date(), "%Y.%m.%d"), ".", "rarefied.phyloseq_melted_wide.txt"), col.names = NA, sep = "\t")


df_melted_rare %>%
  dplyr::rename(ASVID = OTU) %>%
  write.table(., file = paste0(project_name, ".", format(Sys.Date(), "%Y.%m.%d"), ".", "rarefied.phyloseq_melted_long.txt"), col.names = NA, sep = "\t")


# Rarefied table followed by relative abundance
psmelt(final_phyloseq_rarefied_16s %>% transform_sample_counts(.,function(x) (x / sum(x))*100)) %>%
  dplyr::rename(ASVID = OTU) %>%
  write.table(., file = paste0(project_name, ".", format(Sys.Date(), "%Y.%m.%d"), ".", "rarefied.relative_phyloseq_melted_long.txt"), col.names = NA, sep = "\t")


psmelt(final_phyloseq_rarefied_16s %>% transform_sample_counts(.,function(x) (x / sum(x))*100)) %>%
  dplyr::rename(ASVID = OTU) %>%
  tidyr::pivot_wider(names_from = "Sample", values_from = "Abundance") %>%
  write.table(., file = paste0(project_name, ".", format(Sys.Date(), "%Y.%m.%d"), ".", "rarefied.relative_phyloseq_melted_wide.txt"), col.names = NA, sep = "\t")



# Output Alpha diversity rarefied ----
output_alpha_div(otu_table(final_phyloseq_rarefied_16s), project = project_name, is_rarefied = TRUE)


# Output beta diversity rarefied ----
output_beta_div(otu_table(final_phyloseq_rarefied_16s), project = project_name, is_rarefied = TRUE)

# Beta diversity using unifrac distance ----
# Since UniFrac distance has not been implemented in the vegan package, we have to use Phyloseq package to derive it 

write.table(as.matrix(phyloseq::distance(final_phyloseq_rarefied_16s, method = "wunifrac")),
   file = paste0(project_name, ".", format(Sys.Date(), "%Y.%m.%d"), ".","rarefied.weighted_unifrac_dist.txt"), col.names = NA, sep = "\t"
 )



# Log session info ---
sink(paste0(project_name, ".", format(Sys.Date(), "%Y.%m.%d"), ".sessionInfo.txt"))
sessionInfo()
cat("\n \n set.seed is 711")
sink()
