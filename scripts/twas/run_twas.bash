
#This script runs twas analysis - later combined to main rmd

path="scripts/twas"
data="data/twas"
pos_file_names="$path/pos_file_names.txt"
sumstat_dir="sumstats_r12"
out_dir="results_r12"
name="RES_HT"
#name=$1

while read -r line 
do
	#Name of inputfile and its line number and name of tissue extracted
  posfile=$line
	tmp=${posfile#*GTExv8.*.}
	tissue=${tmp%.pos}
	echo $tissue
	  
	#In all chromosomes
	for chr in {1..22}
  do
     Rscript $path/fusion_twas-master/FUSION.assoc_test.R \
				--sumstats $data/$sumstat_dir/${name}_fortwas.gz \
				--weights $data/weights/$posfile \
				--weights_dir $data/weights \
				--ref_ld_chr $data/LDREF/1000G.EUR. \
				--chr $chr \
				--out $data/$out_dir/$name/$name.$tissue.twas.$chr.dat &


  done
  wait
done <$pos_file_names

#Combine and filter results
#62300 for all genes in selected tissues
#Range in selected tissues: 1250 (Kidney-cortex) - 9905 (Adipose)

cat $data/$out_dir/$name/$name.*.twas.*.dat | awk 'NR == 1 || $NF < 0.05/9905' \
 > $data/$out_dir/$name.twas.summary2.tsv


