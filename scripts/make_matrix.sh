#!/bin/bash

num_columns=`cat uniq.txt | wc -l`

declare -A matrix
num_rows=192

for ((i=1;i<=num_columns;i++)); do
  for ((j=1;j<=num_rows;j++)); do
    matrix[$j,$i]=0
  done
done

for i in `seq 1 $num_columns`; do
  line_num=0
  while read line; do
    #echo $line
    stat=`echo $line | awk '{print $11}'`
    snv=`echo $line | awk '{print $5}'`
    name_file=`echo $line | awk '{print $8}' | sed -e 's/;/ /g' | awk '{print $1}'`
    name_comp=`sed ''$i'!d' names.txt`
    if [ $name_file == $name_comp ]; then
      snv_number=0
      while read snv_list; do 
        if [ $snv == $snv_list ]; then
          if [ $stat -le 3 ]; then
            let matrix[$snv_number,$i]=${matrix[$snv_number,$i]}+1
          elif [ $stat -ge 3 ]; then
            let snv_number_tmp=$snv_number+96 
            let matrix[$snv_number_tmp,$i]=${matrix[$snv_number_tmp,$i]}+1
          fi  
        fi
        let snv_number=$snv_number+1
      done < snvs2.txt
    fi
    let line_num=$line_num+1
  done < mutations.txt    
done


f1="%$((${#num_rows}+1))s"
f2=" %9s"

printf "$f1" ''
for ((i=1;i<=num_rows;i++)) do
    printf "$f2" $i
done
echo

for ((j=1;j<=num_columns;j++)) do
    printf "$f1" $j
    for ((i=1;i<=num_rows;i++)) do
        printf "$f2" ${matrix[$i,$j]}
    done
    echo
done
