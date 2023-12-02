
#Downloads each file listed in "namefile"  from the path "url"

namefile=$1
url="https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/EUR"

while read file_name; do
	wget  $url/$file_name &
done <$namefile
