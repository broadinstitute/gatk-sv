#!/bin/bash
image=$1
file=$2

#Describe the Y positions of reads
seq 50 +15 20000 > RANGE


#Extract the start X and Y position of the reads from the file input

StartX=$(cat $file | awk -v OFS="\t" '{if (NR >1) print $1}')
#EndX=$(cat $file | awk -v OFS="\t" '{if (NR >1) print $2}')
#StartX_bottom=$(cat $file | awk -v OFS="\t" '{if (NR >1) print $4}')
#EndX_bottom=$(cat $file | awk -v OFS="\t" '{if (NR >1) print $5}')
Limit=$(cat $file | awk -v OFS="\t" '{if (NR >1) print $NF}')


#Extract the full list of reads
#Y-axis upper haplotype

cat $image  | awk '{ gsub("y=","",$3) }1' | awk '{print $3}' | tr -d '"' | awk -v y="$Limit" '{if ($1>31 && $1<y) print}' | awk '!seen[$0]++' > "$image"_TEMP
cat "$image"_TEMP | while read line;do cat RANGE | grep -w $line; done > "$image"_TEMP2
grep "\S" "$image"_TEMP2 |  awk '!seen[$0]++' > "$image"_TEMP3

#Y-axis bottom haplotype
Limit_bottom=$((Limit+20))
cat $image  | awk '{ gsub("y=","",$3) }1' | awk '{print $3}' | tr -d '"' | awk -v y="$Limit_bottom" '{if ($1>y) print}' | awk '!seen[$0]++' > "$image"_TEMP_bottom
cat "$image"_TEMP_bottom | while read line;do cat RANGE | grep -w $line; done > "$image"_TEMP2_bottom
grep "\S" "$image"_TEMP2_bottom |  awk '!seen[$0]++' > "$image"_TEMP3_bottom


#Upper values
Acount=0
Tcount=0
Gcount=0
Ccount=0

#Bottom values
Acount_bottom=0
Tcount_bottom=0
Gcount_bottom=0
Ccount_bottom=0

#Count nucleotide based on haplotype number

if [ $Limit == 'NA' ] && [ $StartX != 'NA' ]
then
  Lines=$(cat "$image"_TEMP3)
  for line in $Lines
  do
    Yaxis=$(echo $line | awk '{print "y=""\""$1"\""}' | awk '{ gsub (" ", "", $0); print}')
    A=$(cat $image | grep -w $Yaxis | grep -o '>A<' | wc -l)
    T=$(cat $image | grep -w $Yaxis | grep -o '>T<' | wc -l)
    G=$(cat $image | grep -w $Yaxis | grep -o '>G<' | wc -l)
    C=$(cat $image | grep -w $Yaxis | grep -o '>C<' | wc -l)
    Sum=$((A+T+G+C))
    echo $Sum
    Acount=$((Acount+A))
    Tcount=$((Tcount+T))
    Gcount=$((Gcount+G))
    Ccount=$((Ccount+C))
  done >> "$image"_Sum.txt
  Acount_bottom="NA"
  Tcount_bottom="NA"
  Gcount_bottom="NA"
  Ccount_bottom="NA"
  Total=$((Acount+Tcount+Gcount+Ccount))
  median=$(sort -n "$image"_Sum.txt | awk ' { a[i++]=$1; } END { print a[int(i/2)]; }')
  Total_bottom="NA"
  median_bottom="NA"
fi

if [ $Limit != 'NA' ] && [ $StartX != 'NA' ]
then
  Lines=$(cat "$image"_TEMP3)
  for line in $Lines
  do
    Yaxis=$(echo $line | awk '{print "y=""\""$1"\""}' | awk '{ gsub (" ", "", $0); print}')
    A=$(cat $image | grep -w $Yaxis | grep -o '>A<' | wc -l)
    T=$(cat $image | grep -w $Yaxis | grep -o '>T<' | wc -l)
    G=$(cat $image | grep -w $Yaxis | grep -o '>G<' | wc -l)
    C=$(cat $image | grep -w $Yaxis | grep -o '>C<' | wc -l)
    Sum=$((A+T+G+C))
    echo $Sum
    Acount=$((Acount+A))
    Tcount=$((Tcount+T))
    Gcount=$((Gcount+G))
    Ccount=$((Ccount+C))
  done >> "$image"_Sum.txt
  Total=$((Acount+Tcount+Gcount+Ccount))
  median=$(sort -n "$image"_Sum.txt | awk ' { a[i++]=$1; } END { print a[int(i/2)]; }')
fi

if [ $Limit != 'NA' ] && [ $StartX != 'NA' ]
then  
  Lines2=$(cat "$image"_TEMP3_bottom)
  for line2 in $Lines2
  do
    Yaxis_bottom=$(echo $line2 | awk '{print "y=""\""$1"\""}' | awk '{ gsub (" ", "", $0); print}')
    A_bottom=$(cat $image | grep -w $Yaxis_bottom | grep -o '>A<' | wc -l)
    T_bottom=$(cat $image | grep -w $Yaxis_bottom | grep -o '>T<' | wc -l)
    G_bottom=$(cat $image | grep -w $Yaxis_bottom | grep -o '>G<' | wc -l)
    C_bottom=$(cat $image | grep -w $Yaxis_bottom | grep -o '>C<' | wc -l)
    Sum_bottom=$((A_bottom+T_bottom+G_bottom+C_bottom))
    echo $Sum_bottom
    Acount_bottom=$((Acount_bottom+A_bottom))
    Tcount_bottom=$((Tcount_bottom+T_bottom))
    Gcount_bottom=$((Gcount_bottom+G_bottom))
    Ccount_bottom=$((Ccount_bottom+C_bottom))
  done >> "$image"_Sum_bottom.txt
  Total_bottom=$((Acount_bottom+Tcount_bottom+Gcount_bottom+Ccount_bottom))
  median_bottom=$(sort -n "$image"_Sum_bottom.txt | awk ' { a[i++]=$1; } END { print a[int(i/2)]; }')
fi  
  
if [ $Limit == 'NA' ] && [ $StartX == 'NA' ]
then
  Total="NA"
  median="NA"
  Total_bottom="NA"
  median_bottom="NA"
fi

printf '%s\t%s\t%s\t%s\n%s\t%s\t%s\t%s\n' Total_interrupting_nucleotides_upper Median_upper Total_interrupting_nucleotides_bottom Median_bottom $Total $median $Total_bottom $median_bottom

rm -f "$image"_TEMP || true
rm -f "$image"_TEMP2 || true
rm -f "$image"_TEMP3 || true
rm -f "$image"_TEMP_bottom || true
rm -f "$image"_TEMP2_bottom || true
rm -f "$image"_TEMP3_bottom || true
rm RANGE
rm -f "$image"_Sum.txt || true
rm -f "$image"_Sum_bottom.txt || true
