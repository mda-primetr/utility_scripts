

# Requirements 

## Input files
* MetaPhlAn 4 count output: file that ends withs readsStats.table.txt 
* MetaPhlAn 4 relative abundance output: file that contains relative abundance or proportional data
* Project Name : A single word or multiple words separated by "_" (Underscore)
* MetaPhlAn 4 database compatible tree in Newick format 

## R packages
| Package Name | Version | Purpose  | Repository
|:--|:--|:--|:--|
|phyloseq  | >=1.50.0 | Transform and pack microbiome data. Calculate weighted UniFrac distance. | Bioconductor
|readr  | >=2.1.5 | Read BacMap.alpha table |CRAN
|dplyr  | >=1.1.4 | Transform tabular data |CRAN
|tidyr  | >=1.3.1 | Pivot tabular data  |CRAN
|tibble | >=3.2.1 | create tidyverse friendly data frame |CRAN
|magrittr | >=2.0.3 | for %>% pipe |CRAN
|vegan | >=2.6-8 | Alpha diversity indices and beta diversity distances| CRAN
| ape | >=5.8 | Read newick format MetaPhlAn phylogenetic tree|CRAN



# Output files 

#### Phyloseq objects
- Phyloseq object using counts data in the form of RDS file <project_name.YYYY.MM.DD.Metaphlan4.phyloseq_counts.rds>
- Phyloseq object using relative abundance data in the form of RDS file <project_name.YYYY.MM.DD.Metaphlan4.phyloseq_prop.rds>
#### Melted phyloseq objects
- Melted Phyloseq object in the form of wide format tab separated file <project_name.YYYY.MM.DD.Metaphlan4.phyloseq_melted_wide.txt>
- Melted Phyloseq object using proportions or relative abundance in wide format tab separated file <project_name.YYYY.MM.DD.Metaphlan4.phyloseq_prop_melted_wide.txt>
- Melted Phyloseq object in the form of long format tab separated file <project_name.YYYY.MM.DD.Metaphlan4.phyloseq_melted_long.txt>
- Melted Phyloseq object using proportions or relative abundance in long format tab separated file <project_name.YYYY.MM.DD.Metaphlan4.phyloseq_prop_melted_long.txt>
#### Alpha diversity metrics
- Tab delimited file showing number of reads classified into taxonomy and alpha diversity metrics into following columns with rows representing samples <project_name.YYYY.MM.DD.Metaphlan4.no_of_classified_reads_and_alpha_diversity.txt>
- Alpha diversity indices are within sample indices independent of other samples
	- Number of classified reads 
	- Observed number of unique species per sample 
	- InvSimpson showing inverse simpson alpha diversity index
	- Shannon diversity index
#### Beta diversity metrics
- Tab delimited Beta diversity indices or distances (showing relative distances among samples).
	- bray_counts_dist.txt: Pairwise Bray distance based on the counts  <project_name.YYYY.MM.DD.Metaphlan4.bray_counts_dist.txt>
	- bray_relative_dist.txt: Pairwise Bray distance based on the relative abundance data <project_name.YYYY.MM.DD.Metaphlan4.bray_prop_dist.txt>
	- jaccard_dist.txt: Beta diversity distance based on the presence and absence of species regardless of their abundance <project_name.YYYY.MM.DD.Metaphlan4.jaccard_dist.txt>
	- weighted_unifract_dist.txt: Beta diversity distance based on the phylogenetic tree and proportional  <project_name.YYYY.MM.DD.Metaphlan4.prop_weighted_unifrac_dist.txt>


####  R session information
- R session Log <project_name.YYYY.MM.DD.Metaphlan4.sessionInfo.txt>
	- The version of R you're running
	- Your operating system details
	- The base packages that are loaded
	- All other attached packages and their versions
	- Loaded packages that aren't attached

