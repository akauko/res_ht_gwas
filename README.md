# res_ht_gwas
**Code for the project: Genome-Wide Association Study of Apparent Treatment-Resistant Hypertension**

We conducted a genome-wide association analysis for 39737 cases of apparent treatment resistant hypertension and  15996 individuals with mild hypertension. Top varints were replicate at UK Biobank. aTRH-risk loci were characterised using pathway enrichment, a transcriptome-wide association study (TWAS), and Mendelian randomization (MR).

**FinnGen**

* Data: FinnGen https://www.finngen.fi/en
* Description of variables: https://risteys.finngen.fi/
* Finngen Data Specifics: https://finngen.gitbook.io/finngen-handbook/finngen-data-specifics
* Data access: https://site.fingenious.fi/en/
  
**UKBB**

* UKBB Variables: https://biobank.ndph.ox.ac.uk/showcase/index.cgi
* UKBB Research analysis platform: https://www.ukbiobank.ac.uk/enable-your-research/research-analysis-platform

**Methods**

* Regenie: https://rgcgithub.github.io/regenie/
* Susie: https://stephenslab.github.io/susieR/articles/finemapping.html
* VEGAS2: https://vegas2.qimrberghofer.edu.au/
* TWAS/FUSION: http://gusevlab.org/projects/fusion/
* TwoSampleMR: https://mrcieu.github.io/TwoSampleMR/index.html


```
ht_prs_preg
├── README.md		# Overview
├── res_ht_gwas3.rmd 	# R markdown for the analysis
├── res_ht_gwas3.html	# html generated from the R markdown
├── scripts
	├── functions.R	 	# Minor R functions, required by res_ht_gwas3.rmd
	├── gen_res_r12.R	# Creates medication use variables, required by res_ht_gwas3.rmd
	├── fg_pheno_short.txt	# List of used phenotype variables, required by res_ht_gwas3.rmd
	├── rsdify.py 		# Adds rsid:s to FinnGen styled GWAS summaries, provied by FinnGen
	├── README_rsdify.MD 	# Readme for rsdify.MD
	├── vegas2 		# Direct for vegas2 scripts
		├── run_vegas2.bash	# Runs vegas
		├── vegas2v2.pl		# Vegas script (I have added possibility to give directory paths)
		├── generate_GOs.R	# Legacy script for creating own pathway definitons, not used.
	├── twas 		# Directory for twas scripts
		├── run_vegas2.bash	# Runs vegas
		├── vegas2v2.pl		# Vegas script (I have added possibility to give directory paths)
		├── generate_GOs.R	# Legacy script for creating own pathway definitons, not used.
	

├── regenie_input_resht12
	├── blaa		                   #
├── data
	├── blaa		                   #
├── ukbb                         #
	├── blaa
├── local
   ├── blaa		

```
