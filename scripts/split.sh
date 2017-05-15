while read filename; do 
  echo $filename
  grep $filename mutations.txt > split/list_$filename.txt
done < names.txt