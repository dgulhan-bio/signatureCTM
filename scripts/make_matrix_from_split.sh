declare -a array_types=( $(for i in {1..192}; do echo 0; done) )

for file in split/batch$1/*; do
  echo $file
  filename=`echo $file | sed 's/\// /g' | awk '{print $3}'`
  echo 'filename' $filename
  while read line; do
    stat=`echo $line | awk '{print $11}'`
    snv=`echo $line | awk '{print $5}'`
    snv_number=0
    while read snv_list; do
      if [ $snv == $snv_list ]; then
        if [ $stat -le 3 ]; then
          let array_types[$snv_number]=${array_types[$snv_number]}+1
        elif [ $stat -ge 3 ]; then
          let snv_number_tmp=$snv_number+96
          let array_types[$snv_number_tmp]=${array_types[$snv_number_tmp]}+1
        fi
      fi
      let snv_number=$snv_number+1
    done < snvs2.txt
  done < $file
  for i in `seq 0 191`; do
    echo ${array_types[$i]} >> split_out/list_$filename
  done
done



