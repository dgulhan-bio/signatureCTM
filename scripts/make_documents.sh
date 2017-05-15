for file in split/batch$1/*; do
  echo $file
  filename=`echo $file | sed 's/\// /g' | awk '{print $3}'`
  echo 'filename' $filename
  while read line; do
    stat=`echo $line | awk '{print $11}'`
    snv=`echo $line | awk '{print $5}'`
    snv_number=1
    while read snv_list; do
      if [ $snv == $snv_list ]; then
        if [ $stat -le 3 ]; then
          word=`sed ''$snv_number'!d' snv_list_documents.txt`
#          echo $word
          word=$word's'
#          echo $word
          echo $word >> split_out_documents_2/$filename
        elif [ $stat -ge 3 ]; then
          word=`sed ''$snv_number'!d' snv_list_documents.txt`
          word=$word'd'
          echo $word >> split_out_documents_2/$filename
        fi
      fi
      let snv_number=$snv_number+1
    done < snvs2.txt
  done < $file
done
