#!/bin/sh

for file in *
  do awk '{print $6}' $file > snv_$file
done

mkdir snv
mv snv_* snv
rm *
mv snv/* .
rm -r snv

for file in *
  do sed -e 's/\[/p/g' $file| sed -e 's/\]/p/g' | sed -e 's/>/to/g' | awk '{print tolower($0)}' > lower_$file
done

for file in *
  do for i in {1..96}
    do alp=`sed ''$i'!d' ../vocabAlp.txt`
    first=`sed ''$i'!d' ../vocabNoPunc.txt`
    echo $alp $first
    sed -e 's/'$first'/'$alp'/g' $file > tmp_$file
    mv -f tmp_$file $file
 done
done

rm *.sh