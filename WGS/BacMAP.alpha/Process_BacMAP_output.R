#!/usr/bin/env Rscript

# Process_BacMap_output.R
#
# Usage:
#   Rscript Process_BacMap_output /path/to/BacMap_output_table Project_name
#
# Reads a BacMap output table file and constructs a phyloseq object. Saves the phyloseq object
# as an RDS file for later use.
# Rest of the packages should be ok

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    cat("Usage: Rscript Process_BacMap_output <input_tsv> <project name>\n")
    quit(status = 1)
}

bacmap_file <- args[1]
project_name <- args[2]

# Load required packages
suppressPackageStartupMessages({
    library(vegan)
    library(phyloseq)
    library(dplyr)
    library(tidyr)
    library(readr)
    library(tibble)
    library(magrittr) # %>% pipe symbol
})

# Read the raw data ----
d1 <- read_tsv(bacmap_file)

# Get the taxa part
# Format the taxa
d1_tax <- d1 %>%
    dplyr::select(Taxa) %>%
    separate(Taxa, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "TaxID"), sep = ";", convert = TRUE) %>%
    mutate(across(Kingdom:Species, ~ gsub("superkingdom=|phylum=|class=|order=|family=|genus=|species=", "", .x))) %>%
    mutate(TaxID = paste0("TaxID", 1:nrow(d1) - 1)) %>%
    column_to_rownames(var = "TaxID") %>%
    mutate(across(Kingdom:Species, ~ gsub("taxId=[0-9]*", NA, .x))) %>%
    mutate(across(Kingdom:Species, ~ gsub("NA", "", .x))) %>%
    mutate(Kingdom = ifelse(Kingdom == "", "Unknown", Kingdom)) %>%
    mutate(Phylum = ifelse(!is.na(Phylum), paste0("p__", Phylum), Phylum)) %>%
    mutate(Class = ifelse(!is.na(Class), paste0("c__", Class), Class)) %>%
    mutate(Order = ifelse(!is.na(Order), paste0("o__", Order), Order)) %>%
    mutate(Family = ifelse(!is.na(Family), paste0("f__", Family), Family)) %>%
    mutate(Genus = ifelse(!is.na(Genus), paste0("g__", Genus), Genus)) %>%
    mutate(Species = ifelse(!is.na(Species), paste0("s__", Species), Species)) %>%
    rowwise() %>%
    mutate(Phylum = case_when(
        is.na(Phylum) ~ last(c_across(Kingdom)[!is.na(c_across(Kingdom))]),
        TRUE ~ as.character(Phylum)
    )) %>%
    mutate(Class = case_when(
        is.na(Class) ~ last(c_across(Phylum)[!is.na(c_across(Phylum))]),
        TRUE ~ as.character(Class)
    )) %>%
    mutate(Order = case_when(
        is.na(Order) ~ last(c_across(Class)[!is.na(c_across(Class))]),
        TRUE ~ as.character(Order)
    )) %>%
    mutate(Family = case_when(
        is.na(Family) ~ last(c_across(Order)[!is.na(c_across(Order))]),
        TRUE ~ as.character(Family)
    )) %>%
    mutate(Genus = case_when(
        is.na(Genus) ~ last(c_across(Family)[!is.na(c_across(Family))]),
        TRUE ~ as.character(Genus)
    )) %>%
    mutate(Species = case_when(
        is.na(Species) ~ last(c_across(Genus)[!is.na(c_across(Genus))]),
        TRUE ~ as.character(Species)
    )) %>%
    mutate(Phylum = case_when(
        !grepl("p__", Phylum) ~ paste0("LKT_", Phylum),
        TRUE ~ as.character(Phylum)
    )) %>% # This is for labelling unknown taxa levels
    mutate(Class = case_when(
        !grepl("c__", Class) ~ paste0("LKT_", Class),
        TRUE ~ as.character(Class)
    )) %>%
    mutate(Order = case_when(
        !grepl("o__", Order) ~ paste0("LKT_", Order),
        TRUE ~ as.character(Order)
    )) %>%
    mutate(Family = case_when(
        !grepl("f__", Family) ~ paste0("LKT_", Family),
        TRUE ~ as.character(Family)
    )) %>%
    mutate(Genus = case_when(
        !grepl("g__", Genus) ~ paste0("LKT_", Genus),
        TRUE ~ as.character(Genus)
    )) %>%
    mutate(Species = case_when(
        !grepl("s__", Species) ~ paste0("LKT_", Species),
        TRUE ~ as.character(Species)
    )) %>%
    ungroup() %>%
    mutate(TaxID = paste0("TaxID", 1:nrow(d1) - 1)) %>%
    column_to_rownames(var = "TaxID") %>%
    as.matrix()


# Get the OTU count parts
d1_counts <- d1 %>%
    dplyr::select(-Taxa) %>%
    mutate(TaxID = paste0("TaxID", 1:nrow(d1) - 1)) %>%
    column_to_rownames(var = "TaxID") %>%
    as.matrix()

OTU <- otu_table(d1_counts, taxa_are_rows = TRUE)
TAX <- tax_table(d1_tax)

physeq <- phyloseq(OTU, TAX)


# Save the phyloseq object as an RDS file
saveRDS(physeq, file = paste0(project_name, ".", format(Sys.Date(), "%Y.%m.%d"), ".BacMap.phyloseq.rds"))





# function  Alpha diversity ----
output_alpha_div <- function(OTU_mat, project = project_name) {
    write.table(
        data.frame(
            cbind(
                number_of_classified_reads = colSums(OTU_mat),
                Observed = colSums(OTU_mat > 0),
                InvSimpson = vegan::diversity(t(OTU_mat), index = "invsimpson"),
                Shannon = vegan::diversity(t(OTU_mat), index = "shannon")
            )
        ),
        file = paste0(project, ".", format(Sys.Date(), "%Y.%m.%d"), ".BacMap.no_of_classified_reads_and_alpha_diversity.txt"),
        col.names = NA, sep = "\t"
    )
}


# function Beta diversity ----

output_beta_div <- function(OTU_mat, project = project_name) {
    # Bray distance as is
    write.table(as.matrix(vegan::vegdist(t(OTU_mat), method = "bray")),
        file = paste0(project, ".", format(Sys.Date(), "%Y.%m.%d"), ".BacMap.bray_counts_dist.txt"), col.names = NA, sep = "\t"
    )

    # Bray distance on relative abundances
    write.table(as.matrix(vegan::vegdist(apply(t(OTU_mat), 2, function(x) (x / sum(x))), method = "bray")),
        file = paste0(project, ".", format(Sys.Date(), "%Y.%m.%d"), ".BacMap.bray_relative_dist.txt"), col.names = NA, sep = "\t"
    )

    # Jaccard distance
    write.table(as.matrix(vegan::vegdist(t(OTU_mat), method = "jaccard")),
        file = paste0(project, ".", format(Sys.Date(), "%Y.%m.%d"), ".BacMap.jaccard_dist.txt"), col.names = NA, sep = "\t"
    )
}


# Output melted phyloseq object ----
psmelt(physeq) %>%
    tidyr::pivot_wider(names_from = "Sample", values_from = "Abundance") %>%
    dplyr::rename(TaxID = OTU) %>%
    write.table(., file = paste0(project_name, ".", format(Sys.Date(), "%Y.%m.%d"), ".BacMap.phyloseq_melted_wide.txt"), col.names = NA, sep = "\t")



psmelt(physeq) %>%
    dplyr::rename(TaxID = OTU) %>%
    write.table(., file = paste0(project_name, ".", format(Sys.Date(), "%Y.%m.%d"), ".BacMap.phyloseq_melted_long.txt"), col.names = NA, sep = "\t")



# Output Alpha diversity ----
output_alpha_div(d1_counts, project = project_name)

# Output beta diversity ----
output_beta_div(d1_counts, project = project_name)


# Log session info ---
sink(paste0(project_name, ".", format(Sys.Date(), "%Y.%m.%d"), ".BacMap.sessionInfo.txt"))
sessionInfo()
sink()
