README for BacMAP Output:

1. Requirements 
BacMAP count output file: BacMAP output file that ends withs readsStats.table.txt 
Project Name : A single word or multiple words separated by "_" (Underscore)


R packages:
- phyloseq_1.50.0 (Package for analyze and transform microbiome data for further analysis)
- readr_2.1.5 (package to read BacMAP table)
- dplyr_1.1.4 (package to transform data)
- tidyr_1.3.1 (package to pivot data)
- tibble_3.2.1 (package to create dataframe)
- magrittr_2.0.3 (pacakge to use %>% pipe)
- vegan_2.6-8 (Package to calculate alpha diversity indices and beta diversity pairwise distance)

2. Output files 
- Phyloseq object in the form of RDS file <project_name.YYYY.MM.DD.BacMAP.phyloseq.rds>
- Melted Phyloseq object in wide format <project_name.YYYY.MM.DD.BacMAP.phyloseq_melted_wide.txt>
- Melted Phyloseq object in long format <project_name.YYYY.MM.DD.BacMAP.phyloseq_melted_long.txt>
- Tab delimited file showing number of reads classified into taxonomy and alpha diversity metrics into following columns with rows representing samples <project_name.YYYY.MM.DD.BacMAP.no_of_classified_reads_and_alpha_diversity.txt>
- Alpha diversity indices are within sample indices independent of other samples
	- Number of classified reads 
	- Observed number of unique species per sample 
	- InvSimpson showing inverse simpson alpha diversity index
	- Shannon diversity index
- Tab delimited Beta diversity indices or distances (showing relative distances among samples).
- Beta diversity indices or distances are dependent on other samples
	- bray_counts_dist.txt: Pairwise Bray distance based on the counts  <project_name.YYYY.MM.DD.BacMAP.bray_counts_dist.txt>
	- bray_relative_dist.txt: Pairwise Bray distance based on the relative abundance data <project_name.YYYY.MM.DD.BacMAP.bray_relative_dist.txt>
	- jaccard_dist.txt: Beta diversity distance based on the presence and absence of Species regardless of their abundance <project_name.YYYY.MM.DD.BacMAP.jaccard_dist.txt>


3. R session Log <project_name.YYYY.MM.DD.BacMAP.sessionInfo.txt>
- The version of R you're running
- Your operating system details
- The base packages that are loaded
- All other attached packages and their versions
- Loaded packages that aren't attached