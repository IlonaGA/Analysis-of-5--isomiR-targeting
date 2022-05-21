# Analysis-of-5-isomiR-targeting

### Background  
MiRNAs are short non-coding molecules that are able to regulate gene expression post-transcriptionally (1). Previous studies showed that miRNAs may have variabilities at their 5ʹ- and 3ʹ-ends, which appear during miRNA proccessing due to innacurate cleavage of Drosha and Dicer enzymes, these variants of miRNAs were named isomiRs (2, 3). Each miRNA binds to its target-gene via seed region, that is nucleotides sequence located at positions 2-7 of 5'-end. Single nucleotide variation at the 5ʹ-end could alter the seed sequence, which can lead to different target gene binding.
MiRNA processing (4):  
<img src="[https://github.com/IlonaGA/Analysis-of-5--isomiR-targeting/blob/main/Images/res.png](https://www.mdpi.com/ijms/ijms-21-01723/article_deploy/html/images/ijms-21-01723-g001.png)" width="500"/> 


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

**5_plots.py** draws scatterplots TCGA-CHOL, TCGA-PCPG, these two cancers have the biggest and the smallest Spearmanʹs correlation coefficients between the number of predicted targets and the number of anti-correlating targets. Scrips also creates a dotplot which shows Spearmanʹs correlation coefficient for each cancer.  

Plot:  
<img src="https://github.com/IlonaGA/Analysis-of-5--isomiR-targeting/blob/main/Images/res.png" width="800"/>  


### Results:  
- Activity of 5’-isomiRs adjusted for background correlation effects was infered;  
- The 5’-isomiR-gene interaction networks for each cancer from TCGA were constructed;  
- Dependency between the number of predicted targets and the number of anti-correlating targets was examined.

### Conclusions:
- The set of microRNAs and their target-genes pairs differs between cancer types, which may have a functional effect;
- Isoforms of one miRNA can have different target genes;
- Distribution of Spearmanʹs correlation coefficients between the number of predicted targets and the number of anti-correlating targets varies for different types of cancer.

### References

1. Hobert,O. (2008) Gene regulation by transcription factors and microRNAs. Science, 319, 1785–6.
2. Morin,R.D., O’Connor,M.D., Griffith,M., Kuchenbauer,F., Delaney,A., Prabhu,A.L., Zhao,Y., McDonald,H., Zeng,T., Hirst,M., et al. (2008) Application of massively parallel sequencing  to microRNA profiling and discovery in human embryonic stem cells. Genome Res., 18, 610–621.
3. Neilsen,C.T., Goodall,G.J. and Bracken,C.P. (2012) IsomiRs – the overlooked repertoire in the dynamic microRNAome. Trends Genet., 28, 544–549.
4. Ali Syeda, Z., Langden, S., Munkhzul, C., Lee, M., & Song, S. J. (2020). Regulatory Mechanism of MicroRNA Expression in Cancer. International journal of molecular sciences, 21(5), 1723. 
