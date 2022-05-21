# Analysis-of-5-isomiR-targeting
### Aim 
To analyse and identify 5’-isomiRs with significant targeting activity in a given set of samples with available isomiR/mRNA expression data.

### Repository description:
Repository contains scripts which were used for calculating miRNAs activity.   
**1_create_raw_and_summary.py** takes RPM-normalized TCGA data and the list of predicted targets for each TCGA project. The result of script is two tables: TCGA_project_name_raw.tsv and TCGA_project_name_summary.tsv.  
  
Example of TCGA_project_name_raw.tsv table:  
  
| isomiR | gene	| isomiR_expression_median | gene_expression_median |	corr | sample size | proba |  
| ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- |  
| hsa-miR-381-3p/0 | AADAC | 8.90024574856556 |	0.845765094812417 |	-0.14747205 |	79 | 0.08140444164119565 |  
| hsa-miR-182-5p/0	| AAGAB	| 13.013437774635646 |	3.59722347655838 | 0.108154766 | 79 | 0.00013937321980098482 |    
  
Where isomiR is miRNA name, gene is target gene, isomiR_expression_median is the median of miRNA expression, gene_expression_median is target gene expression median, corr is Spearman correlation between miRNA expression and target gene expression, sample size is the amount of samples for each miRNA, proba is cdf of Spearman correlation distribution.

Example of TCGA_project_name_summary.tsv table:  

| isomiR | expression_median | predicted_targets |	sign_corr |  
| ------------- | ------------- | ------------- | ------------- |
| hsa-miR-96-5p/0 |	3.877957395587893 |	182 |	19 |  
hsa-miR-424-5p/0 |	9.549397891501824 |	206 |	22 |  

Where isomiR is miRNA name, expression_median is the median of miRNA expression, predicted_targets is number of predicted targets, sign_corr is amount of anti-correlated target genes.

**2_activity_analysis.py** takes result files obtained by **1_create_raw_and_summary.py** and highly expressed miRNAs list an calculates activity for each miRNA.  

Example of of result file table:  

| isomiR | expression_median | activity |	p_value |	FDR |  
| ------------- | ------------- | ------------- | ------------- | ------------- |    
| hsa-miR-29c-3p/0 | 12.386809776223426 |	7.0 |	7.474732213013041e-09 |	4.0363553950270426e-07 |  
| hsa-miR-30a-5p/0 |	17.094263109900616 |	1.0 |	0.011895843208080804 |	0.3211877666181817 |  
  
Where isomiR is miRNA name, expression_median is the median of miRNA expression, activity is miRNA activity, p_value is p-value for activity, FDR is False Discovery Rate.

**3_target_intersection.py** takes files with activities and list of miRNAs and then finds target intersections for them. Script is useful for understanding which miRNAs have the same targets genes.  

**4_interaction_networks.py** requeres raw and summary files and the list of  highly expressed miRNAs. Draws interaction graphs using Cytoscape.

TCGA-COAD example:   
<img src="https://github.com/IlonaGA/Analysis-of-5--isomiR-targeting/blob/main/Images/TCGA-COAD.png" width="800"/>  


### Results:  
- Activity of 5’-isomiRs adjusted for background correlation effects was infered;  
- The 5’-isomiR-gene interaction networks for each cancer from TCGA were constructed;  
- Dependency between the number of predicted targets and the number of anti-correlating targets was examined.

 
