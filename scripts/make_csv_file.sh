for i in `seq 1 192`; do
  touch tmp.txt
  for file in split_out/*; do
    number=`cat $file | sed ''$i'!d'`
    echo $number >> tmp.txt
#    echo $number
  done
  cat tmp.txt | sed -e ':a' -e 'N' -e '$!ba' -e 's/\n/,/g' >> matrix_from_split.txt
  rm tmp.txt
done