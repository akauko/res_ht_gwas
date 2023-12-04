#!/bin/bash
# NOTE : Quote it else use array to avoid problems #
#dir=$1
dir="../data/risk_traits/*.rds"
for file_in in $dir
do
  file_out="${file_in%.*}.tsv.gz"
	echo "Converting $file_in to $file_out..."
  echo
	R --slave -e "readRDS('$file_in') |> data.table::fwrite('$file_out', sep='\t')" &

done
