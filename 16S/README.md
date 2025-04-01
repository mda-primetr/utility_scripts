These files are generated using PRIME-TR 16S production pipeline in collaboration with CCSG Microbiome Research Core. 

README for 16S Output:

1. Requirements 
Biom file: 16S data biom file containing Taxonomy, OTU/ASV counts, sequences and 

Project Name : A single word or multiple words separated by "_" (Underscore)

R packages:
phyloseq_1.50.0 (To analyze and transform microbiome data for further analysis)
rbiom_1.0.3 (To read biom files and it must be this version)
vegan_2.6-8 (To calculate alpha diversity indices and beta diversity pairwise distance)
Biostrings_2.74.0 (To integrate sequences into phyloseq object
plyr_1.1.4 (package to transform data)
tidyr_1.3.1 (package to pivot data)
magrittr_2.0.3 (pacakge to use %>% pipe)

2. Output files 
- Phyloseq object in the form of RDS file <project_name.YYYY.MM.DD.phyloseq.rds>
- Phyloseq object in the form of wide format tab separated file but does not contain sequences and phylogenetic tree <project_name.YYYY.MM.DD.phyloseq_melted_wide.txt>
- Phyloseq object in the form of long format tab separated file but does not contain sequences and phylogenetic tree <project_name.YYYY.MM.DD.phyloseq_melted_long.txt>
- Tab delimited file showing number of reads classified into taxonomy and alpha diversity metrics into following columns with rows representing samples <project_name.YYYY.MM.DD.no_of_classified_reads_and_alpha_diversity.txt>
- Alpha diversity indices are within sample indices independent of other samples
	- Number of classified reads 
	- Observed number of unique ASVs per sample 
	- InvSimpson showing inverse simpson alpha diversity index
	- Shannon diversity index
- Tab delimited Beta diversity indices or distances (showing relative distances among samples).
- Beta diversity indices or distances are dependent on other samples
	- bray_counts_dist.txt: Pairwise Bray distance based on the counts  <project_name.YYYY.MM.DD.bray_counts_dist.txt>
	- bray_relative_dist.txt: Pairwise Bray distance based on the relative abundance data <project_name.YYYY.MM.DD.bray_relative_dist.txt>
	- jaccard_dist.txt: Beta diversity distance based on the presence and absence of ASVs regardless of their abundance <project_name.YYYY.MM.DD.jaccard_dist.txt>

3. Rarefied data
Rarefied microbiome data is a normalization technique which by default reduces the number of taxonomy classified reads in all samples equal to the sample with minimum number of reads.  
This randomly selects the ASVs from the pre-rarefied data so please note the random seed set in the script for reproducibilty.
All rarefied output files have "rarefied" keyword in them and the number of rarefied output files should be same as without rarefaction.


4. R session Log <project_name.YYYY.MM.DD.sessionInfo.txt>
- The version of R you're running
- Your operating system details
- The base packages that are loaded
- All other attached packages and their versions
- Loaded packages that aren't attached
- Random seed that was set for rarefaction function 

