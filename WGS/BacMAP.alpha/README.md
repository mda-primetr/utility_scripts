

# Requirements

## Input files
* BacMAP count output file: BacMAP output file that ends withs readsStats.table.txt 
* Project Name : A single word or multiple words separated by "_" (Underscore)
## R packages

| Package Name | Version | Purpose  | Repository
|:--|:--|:--|:--|
|phyloseq  | >=1.50.0 | Transform and pack microbiome data | Bioconductor
|readr  | >=2.1.5 | Read BacMap.alpha table | CRAN
|dplyr  | >=1.1.4 | Transform tabular data | CRAN
|tidyr  | >=1.3.1 | Pivot tabular data  | CRAN
| tibble | >=3.2.1 | create tidyverse friendly data frame | CRAN
| magrittr | >=2.0.3 | for %>% pipe | CRAN
| vegan | >=2.6-8 | Alpha diversity indices and beta diversity distances| CRAN

***Use R >=4.0***

# Output 

#### Phyloseq object
- Phyloseq object in the form of RDS file <project_name.YYYY.MM.DD.BacMAP.phyloseq.rds>
#### Melted phyloseq objects
- Using counts
	- Melted Phyloseq object in wide format <project_name.YYYY.MM.DD.BacMAP.phyloseq_melted_wide.txt>
	- Melted Phyloseq object in long format <project_name.YYYY.MM.DD.BacMAP.phyloseq_melted_long.txt>
- Using relative proportion
	- Melted Phyloseq object in wide format <project_name.YYYY.MM.DD.BacMap.phyloseq_relative_wide.txt>
	- Melted Phyloseq object in wide format <project_name.YYYY.MM.DD.BacMap.phyloseq_relative_long.txt>
#### Alpha diversity metrics
- Tab delimited file showing number of reads classified into taxonomy and alpha diversity metrics into following columns with rows representing samples <project_name.YYYY.MM.DD.BacMAP.no_of_classified_reads_and_alpha_diversity.txt>
** Alpha diversity indices are within sample indices independent of other samples
	- Number of taxa classified reads 
	- Observed number of unique species per sample 
	- InvSimpson showing inverse simpson alpha diversity index
	- Shannon diversity index
#### Beta diversity metrics
-  Beta diversity indices or distances are dependent on other samples
	- bray_counts_dist.txt: Pairwise Bray distance based on the counts  <project_name.YYYY.MM.DD.BacMAP.bray_counts_dist.txt>
	- bray_relative_dist.txt: Pairwise Bray distance based on the relative abundance data <project_name.YYYY.MM.DD.BacMAP.bray_relative_dist.txt>
	- jaccard_dist.txt: Beta diversity distance based on the presence and absence of Species regardless of their abundance <project_name.YYYY.MM.DD.BacMAP.jaccard_dist.txt>


#### R session information
 - R session Log <project_name.YYYY.MM.DD.BacMAP.sessionInfo.txt>
	- The version of R you're running
	- Your operating system details
	- The base packages that are loaded
	- All other attached packages and their versions
	- Loaded packages that aren't attached