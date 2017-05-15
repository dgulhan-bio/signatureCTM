rm bad.txt
while read line; do
  name=`echo $line | awk '{print $2}'`
  qual=`echo $line | awk '{print $3}'`
#  echo $qual
  if [ ${qual%.*} -ge 1 ]; then
    echo $name >> bad.txt
  fi
done < ALL_SAMPLES.mapd.txt

for name in `cat names.txt | awk '{print $1}' | sed -e 's/\-/ /g' | awk '{print $1}' | uniq`; do
  grep $name bad.txt >> bad_in_this.txt
done
