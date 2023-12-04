
#This script runs twas analysis - later combined to main rmd

#"bash get_qtl file_names.txt" was run at local computer: it fetches all qtl weights
#listed in the 'namefile'. All resulting files were transfered to FinnGen

path="scripts/twas"
data="data/twas"

#All weights in given directeory are unpacked
bash $path/unpack_all_targz.bash $data/weights

#Extract list of snips in LDREF set
#Unpacks all tar.gz files in directory, which is given as parameter

snplist_out="$data/ldref_snplist.txt"
echo "" > $snplist_out
for chr in {1..22}
do
  awk '{print $2}' data/twas/LDREF/*.$chr.bim >> $snplist_out
done

#Number of genes used for analysis from input line numbers  (hopefully)
find $data -name "*EUR*.pos" ! -name "*nofilter*" -exec wc -l {} \; > $data/nr_lines_pos.csv

