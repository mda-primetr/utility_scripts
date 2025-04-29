
# Requirements

## Input files
* Biom file: 16S data biom file containing Taxonomy, OTU/ASV counts, sequences.

* Project Name : A single word or multiple words separated by "_" (Underscore)

## R packages:
| Package Name | Version | Purpose  | Repository
|:--|:--|:--|:--|
|phyloseq  | >=1.50.0 | Transform and pack microbiome data. Calculate weighted UniFrac distance. | Bioconductor
|rbiom | ==1.0.3 | Read biom file | CRAN|
|Biostrings | >=2.74.0 | integrate sequences into phyloseq object | Bioconductor
|readr  | >=2.1.5 | Read BacMap.alpha table |CRAN
|dplyr  | >=1.1.4 | Transform tabular data |CRAN
|tidyr  | >=1.3.1 | Pivot tabular data  |CRAN
|tibble | >=3.2.1 | create tidyverse friendly data frame |CRAN
|magrittr | >=2.0.3 | for %>% pipe |CRAN
|vegan | >=2.6-8 | Alpha diversity indices and beta diversity distances| CRAN


# Output files 

#### Phyloseq objects
- Phyloseq object in the form of RDS file <project_name.YYYY.MM.DD.phyloseq.rds>
- Phyloseq object in the form of RDS file <project_name.YYYY.MM.DD.phyloseq.rarefied.rds>
#### Melted phyloseq objects
- Phyloseq object in the form of wide format tab separated file but does not contain sequences and phylogenetic tree <project_name.YYYY.MM.DD.phyloseq_melted_wide.txt>
- Phyloseq object in the form of long format tab separated file but does not contain sequences and phylogenetic tree <project_name.YYYY.MM.DD.phyloseq_melted_long.txt>
- Phyloseq object in the form of long format tab separated file but does not contain sequences and phylogenetic tree with relative abundance <project_name.YYYY.MM.DD.relative_phyloseq_melted_long.txt>
- Phyloseq object in the form of wide format tab separated file but does not contain sequences and phylogenetic tree with relative abundance <project_name.YYYY.MM.DD.relative_phyloseq_melted_wide.txt>
- Tab delimited file showing number of reads classified into taxonomy and alpha diversity metrics into following columns with rows representing samples <project_name.YYYY.MM.DD.no_of_classified_reads_and_alpha_diversity.txt>
#### Alpha diversity metrics
- Alpha diversity indices are within sample indices independent of other samples
	- Number of classified reads 
	- Observed number of unique ASVs per sample 
	- InvSimpson showing inverse simpson alpha diversity index
	- Shannon diversity index
#### Beta diversity metrics
 Tab delimited Beta diversity indices or distances (showing relative distances among samples).
- Beta diversity indices or distances are dependent on other samples
	- bray_counts_dist.txt: Pairwise Bray distance based on the counts  <project_name.YYYY.MM.DD.bray_counts_dist.txt>
	- bray_relative_dist.txt: Pairwise Bray distance based on the relative abundance data <project_name.YYYY.MM.DD.bray_relative_dist.txt>
	- jaccard_dist.txt: Beta diversity distance based on the presence and absence of ASVs regardless of their abundance <project_name.YYYY.MM.DD.jaccard_dist.txt>
	- weighted_unifrac_dist.txt: Beta diversity distance based on the UniFrac distance <project_name..YYYY.MM.DD.rarefied.weighted_unifrac_dist.txt>

####  Rarefied data
- Rarefied microbiome data is a normalization technique which by default reduces the number of taxonomy classified reads in all samples equal to the sample with minimum number of reads. This randomly selects the ASVs from the pre-rarefied data so please note the random seed set in the script for reproducibilty.
- All rarefied output files have **rarefied** keyword in them and the number of rarefied output files should be same as without rarefaction.

#### R session information
-  R session Log <project_name.YYYY.MM.DD.sessionInfo.txt>
	- The version of R you're running
	- Your operating system details
	- The base packages that are loaded
	- All other attached packages and their versions
	- Loaded packages that aren't attached
	- Random seed that was set for rarefaction function 

