# Analysis-of-5-isomiR-targeting
### Aim 
To analyse and identify 5’-isomiRs with significant targeting activity in a given set of samples with available isomiR/mRNA expression data.

### Repository description:
Repository contains scripts which were used for calculating miRNAs activity. 
- 1_create_raw_and_summary.py takes RPM-normalized TCGA data and the list of predicted targets for each TCGA project. The result of script is two tables: TCGA_project_name_raw.tsv and TCGA_project_name_summary.tsv.  
Example of TCGA_project_name_raw.tsv table:  
  
| isomiR | gene	| isomiR_expression_median | gene_expression_median |	corr | sample size | proba |  
| ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- |  
| hsa-miR-381-3p|0 | AADAC | 8.90024574856556 |	0.845765094812417 |	-0.14747205 |	79 | 0.08140444164119565 |  
| hsa-miR-182-5p|0	| AAGAB	| 13.013437774635646 |	3.59722347655838 | 0.108154766 | 79 | 0.00013937321980098482 |    
  
Where isomiR is miRNA name, gene is target gene, isomiR_expression_median is median for miRNA expression, gene_expression_median is target gene expression median, corr is Spearman correlation between miRNA expression and target gene expression, sample size is the amount of samples for each miRNA, proba is
- 


### Results:
- 38% of highly expressed isomiRs (63 entries) were expressed in 10 or more cancer types, including 11 non-canonical 5’-isomiRs;  
- Was shown that high expression of a non-canonical 5’-isomiR is tied with the high expression of its corresponding canonical form;  
- Activity of 5’-isomiRs adjusted for background correlation effects was infered;  
- The 5’-isomiR-gene interaction networks for each cancer from TCGA were constructed;  
- For a given 5’-isomiR and cancer type broad 5’-isomiR targeting activity was denoted.  

| First Header  | Second Header |
| ------------- | ------------- |
| Content Cell  | Content Cell  |
| Content Cell  | Content Cell  |
