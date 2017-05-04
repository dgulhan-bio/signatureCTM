#!/bin/sh

#strip the list of snvs
for file in *
  do awk '{print $6}' $file > snv_$file
done

mkdir snv
mv snv_* snv
rm *
mv snv/* .
rm -r snv

#remove paranthesis make lower case
for file in *
  do sed -e 's/\[/p/g' $file| sed -e 's/\]/p/g' | sed -e 's/>/to/g' | awk '{print tolower($0)}' > lower_$file
done

#change vocabulary to a simpler format
for file in *
  do for i in {1..96}
    do alp=`sed ''$i'!d' ../vocabAlp.txt`
    first=`sed ''$i'!d' ../vocabNoPunc.txt`
    echo $alp $first
    sed -e 's/'$first'/'$alp'/g' $file > tmp_$file
    mv -f tmp_$file $file
 done
done

rm -f snv_*
rm -f *convert.sh

#remove the first line
for file in *
  do  tail -n +2 $file > test.tmp
  mv -f test.tmp $file
done

touch writeToR.R
echo 'frame3BaseAging<-data.frame(c("' > writeToR.R

nfile=`ls lower_snv_* | wc -l`
count=1
for file in lower*
  do sed -e ':a' -e 'N' -e '$!ba' -e 's/\n/ /g' $file >> writeToR.R
  if [ $count -lt $nfile ]; then
     echo '","' >> writeToR.R
  fi
  let count=$count+1
  echo $count
done
echo '"))' >> writeToR.R

rm -f test.tmp
sed -e ':a' -e 'N' -e '$!ba' -e 's/\n/ /g' writeToR.R > test.tmp
mv -f test.tmp  writeToR.R

echo 'save(frame3BaseAging,file="frame3BaseAging.Rda")' >> writeToR.R