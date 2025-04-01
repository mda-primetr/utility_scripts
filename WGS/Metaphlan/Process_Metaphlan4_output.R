#!/usr/bin/env Rscript

# Process_Metaphlan4_output.R
#
# Usage:
#   Rscript Process_Metaphlan4_output /path/to/metaphlan.readStats.table Project_name Metaphlan4Tree_file
#
# Reads a Metaphlan4 output table file and constructs a phyloseq object. Saves the phyloseq object
# as an RDS file for later use.
# Rest of the packages should be ok
# Metaphlan4Tree is mpa_vJun23_CHOCOPhlAnSGB_202307.nwk for Jun23 database. Please check it.

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    cat("Usage: Rscript Process_Metaphlan4_output <input_tsv> <project name> <tree_file>\n")
    quit(status = 1)
}

metapahlan_file <- args[1] # Metaphlan.read.stats or counts table and *NOT* relative abundance output table
project_name <- args[2]
metaphlan_tree <- args[3]

# Load required packages
suppressPackageStartupMessages({
    library(vegan)
    library(ape)
    library(phyloseq)
    library(dplyr)
    library(tidyr)
    library(readr)
    library(tibble)
    library(magrittr) # %>% pipe symbol
})



# Read the new metaphlan 4 data
t4 <- read_tsv(metapahlan_file) %>%
    separate(col = c(clade), into = c(
        "Kingdom",
        "Phylum",
        "Class",
        "Order",
        "Family",
        "Genus",
        "Species",
        "ID"
    ), sep = "\\|") %>%
    mutate(ID = gsub("t__SGB", "", ID)) %>%
    mutate(ID = gsub("_group", "", ID)) %>%
    filter(!is.na(ID))


t4_tax <- t4 %>%
    dplyr::select(Kingdom:ID) %>%
    dplyr::mutate(Kingdom = ifelse(Kingdom == "", "Unknown", Kingdom)) %>%
    rowwise() %>%
    dplyr::mutate(Phylum = case_when(
        grepl("p__PGB", Phylum) ~ last(c_across(Kingdom)[!is.na(c_across(Kingdom))]),
        TRUE ~ as.character(Phylum)
    )) %>%
    dplyr::mutate(Class = case_when(
        grepl("c__CFGB", Class) ~ last(c_across(Phylum)[!is.na(c_across(Phylum))]),
        TRUE ~ as.character(Class)
    )) %>%
    dplyr::mutate(Order = case_when(
        grepl("o__OFGB", Order) ~ last(c_across(Class)[!is.na(c_across(Class))]),
        TRUE ~ as.character(Order)
    )) %>%
    dplyr::mutate(Family = case_when(
        grepl("f__FGB", Family) ~ last(c_across(Order)[!is.na(c_across(Order))]),
        TRUE ~ as.character(Family)
    )) %>%
    dplyr::mutate(Genus = case_when(
        grepl("g__GGB", Genus) ~ last(c_across(Family)[!is.na(c_across(Family))]),
        TRUE ~ as.character(Genus)
    )) %>%
    dplyr::mutate(Species = case_when(
        grepl("s__GGB", Species) ~ last(c_across(Genus)[!is.na(c_across(Genus))]),
        TRUE ~ as.character(Species)
    )) %>%
    dplyr::mutate(Phylum = case_when(
        !grepl("p__", Phylum) ~ paste0("LKT_", Phylum),
        TRUE ~ as.character(Phylum)
    )) %>% # This is for labelling unknown taxa levels
    dplyr::mutate(Class = case_when(
        !grepl("c__", Class) ~ paste0("LKT_", Class),
        TRUE ~ as.character(Class)
    )) %>%
    dplyr::mutate(Order = case_when(
        !grepl("o__", Order) ~ paste0("LKT_", Order),
        TRUE ~ as.character(Order)
    )) %>%
    dplyr::mutate(Family = case_when(
        !grepl("f__", Family) ~ paste0("LKT_", Family),
        TRUE ~ as.character(Family)
    )) %>%
    dplyr::mutate(Genus = case_when(
        !grepl("g__", Genus) ~ paste0("LKT_", Genus),
        TRUE ~ as.character(Genus)
    )) %>%
    dplyr::mutate(Species = case_when(
        !grepl("s__", Species) ~ paste0("LKT_", Species),
        TRUE ~ as.character(Species)
    )) %>%
    ungroup() %>%
    mutate(SGB = ifelse(grepl("_SGB", Species), Species, paste0(Species, "_SGB", ID))) %>%
    column_to_rownames(var = "ID") %>%
    as.matrix()



t4_otu <- t4 %>%
    dplyr::select(-c(Kingdom:Species)) %>%
    # dplyr::mutate(TaxID = paste0("TaxID", 1:nrow(t4) - 1)) %>%
    column_to_rownames(var = "ID") %>%
    as.matrix()

# Get tree
# In order to merge the tree with the phyloseq object, the tree must have the same tip labels as the rownames of the OTU table
# check if this is the same as the rownames of the OTU table
# b1_tree$tip.label
b1_tree <- ape::read.tree(metaphlan_tree)



# Create phyloseq components ----
OTU <- otu_table(t4_otu, taxa_are_rows = TRUE)
TAX <- tax_table(t4_tax)

physeq <- phyloseq(OTU, TAX, b1_tree)


# Save the phyloseq object as an RDS file
saveRDS(physeq, file = paste0(project_name, ".", format(Sys.Date(), "%Y.%m.%d"), ".Metaphlan4.phyloseq.rds"))



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
        file = paste0(project_name, ".", format(Sys.Date(), "%Y.%m.%d"), ".Metaphlan4.no_of_classified_reads_and_alpha_diversity.txt"),
        col.names = NA, sep = "\t"
    )
}


# function Beta diversity ----

output_beta_div <- function(OTU_mat, project = project_name) {
    # Bray distance as is
    write.table(as.matrix(vegan::vegdist(t(OTU_mat), method = "bray")),
        file = paste0(project_name, ".", format(Sys.Date(), "%Y.%m.%d"), ".Metaphlan4.bray_counts_dist.txt"), col.names = NA, sep = "\t"
    )

    # Bray distance on relative abundances
    write.table(as.matrix(vegan::vegdist(apply(t(OTU_mat), 2, function(x) (x / sum(x))), method = "bray")),
        file = paste0(project_name, ".", format(Sys.Date(), "%Y.%m.%d"), ".Metaphlan4.bray_relative_dist.txt"), col.names = NA, sep = "\t"
    )

    # Jaccard distance
    write.table(as.matrix(vegan::vegdist(t(OTU_mat), method = "jaccard")),
        file = paste0(project_name, ".", format(Sys.Date(), "%Y.%m.%d"), ".Metaphlan4.jaccard_dist.txt"), col.names = NA, sep = "\t"
    )
}


# Output melted phyloseq object ----
psmelt(physeq) %>%
    tidyr::pivot_wider(names_from = "Sample", values_from = "Abundance") %>%
    dplyr::rename(TaxID = OTU) %>%
    write.table(., file = paste0(project_name, ".", format(Sys.Date(), "%Y.%m.%d"), ".Metaphlan4.phyloseq_melted_wide.txt"), col.names = NA, sep = "\t")


psmelt(physeq) %>%
    dplyr::rename(TaxID = OTU) %>%
    write.table(., file = paste0(project_name, ".", format(Sys.Date(), "%Y.%m.%d"), ".Metaphlan4.phyloseq_melted_long.txt"), col.names = NA, sep = "\t")



# Output Alpha diversity ----
output_alpha_div(t4_otu, project = project_name)


# Output beta diversity ----
output_beta_div(t4_otu, project = project_name)


# Log session info ---
sink(paste0(project_name, ".", format(Sys.Date(), "%Y.%m.%d"), ".Metaphlan4.sessionInfo.txt"))
sessionInfo()
sink()
