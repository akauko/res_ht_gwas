#!/bin/bash

#Runs VEGAS2 (https://vegas2.qimrberghofer.edu.au/)
#This is comptationally intensive!  Would be better to do wdl script...

if [ $# -ne 1 ]; then
        echo "please specify name of the variable"
		exit 1
fi

#NAME=$1
#NAME=HYPTENS
NAME=RES_HT
dir_in=data/regenie_r12
dir=data/vegas2
vegas_path=scripts/vegas2
snps_in=$dir_in/${NAME}_rs.gz

#Inputs downloaded from VEGAS2 website
genelocat_in=$dir/glist-hg19.txt
#pathway_ref=$dir/biosystems20160324.vegas2pathSYM
#pathway_ref=$dir/golevel3_biosystemsnoGO.txt
#pathway_ref=$dir/go_level3_vegas.txt
plink_in=/home/ivm/res_ht/data/vegas2/g1000p3_EUR  #VEGAS requires full plink path

#Directories for chromosome-wise files
snps_dir=$dir/snps_in
glists_dir=$dir/glists
genebased_out_dir=$dir/genebased_out

#------------------------------------------------------------------------------

#Gene file and input files are divided by chromosomes and input files are formatted.
for CHR in {1..23} 
do
  grep -E "^$CHR\\s" $genelocat_in >  $glists_dir/glist-hg19.$CHR.txt
  zcat $snps_in | grep -E "^$CHR\\s" | awk 'BEGIN{OFS="\t"}{print $(NF),$5}' | grep rs > $snps_dir/${NAME}_rs_p.$CHR.tsv & 
done
wait

#-----------------------------------------------------------------------------

#VEGAS2 gene-based run - two days at the 'rather big machine' in parallel...
#(bed file does not have X chromosome, so that is not included here)
# wdl would bee more convenient.
for CHR in {1..22} 
do

  perl $vegas_path/vegas2v2.pl -G -snpandp $snps_dir/${NAME}_rs_p.$CHR.tsv \
                               -custom $plink_in \
                               -glist $glists_dir/glist-hg19.$CHR.txt \
                               -out ${NAME}_genebased_$CHR &
  sleep 1m     #To avoid problems with accessing plink file
done

wait
mv ${NAME}_genebased_*.out $genebased_out_dir/

#Combine
head -1 $genebased_out_dir/${NAME}_genebased_1.out > $dir/${NAME}_genebased.out
for CHR in {1..22} 
do
  tail -n+2 $genebased_out_dir/${NAME}_genebased_$CHR.out >> $dir/${NAME}_genebased.out
done

#-----------------------------------------------------------------------------

#VEGAS2 pathway-based run - about one day at the 'rather big machine' for RES_HT and 
#couple of days for HYPTENS

awk '{print $2,$8}' $dir/${NAME}_genebased.out | grep -v Gene|sed 's/"//g'> $dir/${NAME}.geneandp

pathway_ref=$dir/biosystems20160324.vegas2pathSYM
#pathway_ref=$dir/go_level2_vegas.txt
#pathway_ref=$dir/go_level3_vegas.txt

perl $vegas_path/vegas2v2.pl -P -geneandp $dir/${NAME}.geneandp \
                                -geneandpath $pathway_ref \
																-glist $genelocat_in \
																-out ${NAME}_pathway &
															  
wait
#Annoyingly parameter -out didn't print to given directory
mv ${NAME}_pathway.out $dir/

#-------------------------------------------------------------------------------´

#Filter significant genes and pathways









