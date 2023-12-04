# res_ht_gwas

**Code for the project: Genome-Wide Association Study of Apparent Treatment-Resistant Hypertension**

We conducted a genome-wide association analysis apparent treatment resistant hypertension (aRH) and number of medication classes. aRH GWAS with mild hypertension as control was replicated at UK Biobank was selected for downstream analyses (pathway enrichment, a transcriptome-wide association study, and Mendelian randomization).

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
  ├── README.md                         # Overview
  ├── res_ht_gwas3.rmd                  # R markdown for the analysis
  ├── res_ht_gwas3.html                 # html generated from the R markdown
  ├── scripts
      ├── functions.R                       # Minor R functions, required by res_ht_gwas3.rmd
      ├── gen_res_r12.R                     # Creates medication use variables, required by res_ht_gwas3.rmd
      ├── fg_pheno_short.txt                # List of used phenotype variables, required by res_ht_gwas3.rmd
      ├── rsdify.py                         # Adds rsid:s to FinnGen styled GWAS summaries, provied by FinnGen
      ├── README_rsdify.MD                  # Readme for rsdify.MD
      ├── vegas2                                # Direct for vegas2 scripts
          ├── run_vegas2.bash                   # Vegas wrapper
          ├── vegas2v2.pl                       # Vegas script: Modified to enable use of directory paths
      ├── twas                              # Directory for twas scripts
          ├── run_twas.bash                     # Runs twas for given set tissues - parallelized by chromosome
          ├── pos_file_names.txt                # Position file names - specifies tissues of interest, required by run_twas.bash        
          ├── preprocess_for_twas.bash          # Preprocesses files for twas
          ├── unpack_all_targz.bash             # Required by preprocess_for_twas.bash
  ├── regenie_input_resht12                 # Regenie and finemap input (this directory in FinnGen includes also input data)
      ├── RegenieDF12-1.wdl                     # Regenie wdl script, provided by FinnGen
      ├── subwdls.zip                           # Regenie wdl subscripts, provided by FinnGen
      ├── regenie_resht_bin_R12.json            # Regenie parameters for binary variables
      ├── regenie_resht_cont_R12.json           # Regenie parameters for continous variables
      ├── phenolist_resht_bin_regenie.txt       # list of binary variables
      ├── phenolist_resht_cont_regenie.txt      # list of continouos variables
      ├── finemap.wdl                           # Finemapping wdl script, provided by FinnGen
      ├── finemap_sub.wdl.zip                   # Finemapping wdl subscript, provided by FinnGen
      ├── finemap_resht_r12.json                # Finemapping parameters
      ├── phenolist_finemap.txt                 # list of all outcome variables
  ├── data
      ├── ATC-codes_r8_eng_all.csv              # ATC codes present in R8, also problematic codes included; probably old codes
      ├── ATC-codes_final_r8.csv                # ATC codes present in R8, only codes from current Fimea classification included
      ├── ATC-codes_r12.csv                     # ATC codes present in R12, only codes from current Fimea classification included
      ├── <endpoint>_r12.csv                    # OR table for susie hits, calculated from SUSIE summaries and HYPTENS regenie results
      ├── regenie_r12                                                
          ├── <endpoint>_pval_manhattan.png     # Output of regenie pipeline, manhattan plot
          ├── <endpoint>_pval_manhattan_loglog.png   
          ├── <endpoint>_pval_qqplot.png        # Output of regenie pipeline, qqplot plot  
          ├── <endpoint>_summary.txt            # Output of regenie pipeline, top hits
          ├── <endpoint>.SUSIE.cred.summary.tsv # Output of finemapping pipeline, SUSIE summary
      ├── vegas2                            # vegas2 results for RES_HT and HYPTENS
          ├── <endpoint>_genebased.out          # Results for a gene based run
          ├── <endpoint>_genebased_summary.txt  # Significant results for a gene based run
          ├── <endpoint>_pathway.out            # Results for a pathway based run
          ├── <endpoint>_pathway_summary.txt    # Significant results for a pathway based run
      ├── twas
          ├── hgnc_gene_names.txt               # Gene name mapping  
          ├── resulsts_r12                      
              ├── RES_HT.twas.summary2.tsv      # Summary table for twas run, gene names mapped to results  
              ├── resht_twas_r12.jpg            # Plot for twas results
              ├── resht_twas_r12.pdf            # Plot for twas results
  ├── ukbb                          # Scripts and resuts from UKBB RAP 
      ├── create_ukbb_pheno_resht.Rmd   # 1. Create phenotypes for ukbb
      ├── create_ukbb_pheno_resht.html  
      ├── run_liftover.ipynb            # 2. Liftover chip data from GCRh37 to GCRh38, run at ttyd app at UKB RAP
      ├── run_regenie_resht.ipynb       # 3. Run regenie. Includes pre and post processing
      ├── scripts
          ├── draw_manhattan.R          # Required by run_regenie_resht.ipynb
      ├── liftover_plink_beds       # Input for liftover pipeline
          ├── liftover_plink_beds.wdl   # liftover wdl script, provided by UKBB
          ├── liftover_input.json       # liftover parameters
          ├── b37ToHg38.over.chain      # chain file
      ├── data
          ├── ATC-codes_ukbb.csv        # ATC codes used in analysis 
          ├── regenie                   # UKBB regenie results
              ├── manhattan.resht_<variable>.png
              ├── qqplot.resht_<variable>.png            
  ├── local
      ├── scripts                        # Miscellanous scripts run locally 
          ├── run_MR_r12.R                  # Run mendelian randomization for RES_HT and HYPTENS
          ├── risk_traits_plot.R            # Create risk traits plot for RES_HT replicated hits
          ├── check_signif_snip.R           # Check P valus and betas at UKBB data for SUSIE hits    
          ├── filter_ukbb_by_FinnGen_P.R    # Prefiltering step for check_signif_snip.R
          ├── check_ukbb_atc.R              # Checks semimanually created ATC list agains UKBB codes
          ├── get_qtl.bash                  # Fetches eQTL from gtex for TWAS run
          ├── file_names.txt               # List of eQTL file names for TWAS run    
          ├── generate_GOs.R                # Legacy script for creating own pathway definitons for Vegas2, not used!
      ├── data
          ├── atc_all_matches_c0.csv         # C0* codes from https://github.com/PhilAppleby/ukbb-srmed/blob/master/data/atc_all_matches.csv
          ├── ukbb20003_n.csv                # Counts for drug codes https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=20003
          ├── ATC-codes_ukbb.csv             # Semimanually created based on above lists
          ├── replic
              ├── replic_<variabe>_all.csv       # Replication summary
              ├── replic_ukb_combined.csv        # Manually  created: Replicated variants from variables ukbvar and fgvar
      ├── figs
          ├── resht_mr_r12.jpg                # Mendelian randomization plot, RES_HT
          ├── resht_mr_r12.pdf
          ├── hyptens_mr_r12.jpg              # Mendelian randomization plot, HYPTENS
          ├── hyptens_mr_r12.pdf
          ├── resht_risk_r12.jpg              # Risk plot, RES_HT
          ├── resht_risk_r12.pdf    
          ├── manhattan_resht_r12.jpg         # Manhattan with manually added SNP names
          ├── manhattan_resht_r12.afdesign

```
