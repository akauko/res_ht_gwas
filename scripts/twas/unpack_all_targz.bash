#Unpacks all tar.gz files in directory, which is given as parameter

dir=$1

for f in $dir/*.tar.gz
do
	echo $f
  tar -xzvf $f -C $dir &
done


