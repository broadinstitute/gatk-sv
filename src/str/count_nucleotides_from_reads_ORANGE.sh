#!/bin/bash
image=$1
file=$2

seq 15 +10 20000 > RANGE_x
seq 50 +15 20000 > RANGE

StartX=$(cat $file | awk -v OFS="\t" '{if (NR >1) print $1}')
EndX=$(cat $file | awk -v OFS="\t" '{if (NR >1) print $2}')
StartX_bottom=$(cat $file | awk -v OFS="\t" '{if (NR >1) print $4}')
EndX_bottom=$(cat $file | awk -v OFS="\t" '{if (NR >1) print $5}')
Limit=$(cat $file | awk -v OFS="\t" '{if (NR >1) print $NF}')


#Extract the full list of reads
#Y-axis upper haplotype

cat $image  | awk '{ gsub("y=","",$3) }1' | awk '{print $3}' | tr -d '"' | awk -v y="$Limit" '{if ($1>31 && $1<y) print}' | awk '!seen[$0]++' > "$image"_TEMP
cat "$image"_TEMP | while read line;do cat RANGE | grep -w $line; done > "$image"_TEMP2
grep "\S" "$image"_TEMP2 |  awk '!seen[$0]++' > "$image"_TEMP3
#cat RANGE_x | awk -v x="$StartX" -v xe="$EndX" '{if ($1>x && $1<xe) print}' > "$image"_TEMP_X


#Y-axis bottom haplotype
Limit_bottom=$((Limit+20))
cat $image  | awk '{ gsub("y=","",$3) }1' | awk '{print $3}' | tr -d '"' | awk -v y="$Limit_bottom" '{if ($1>y) print}' | awk '!seen[$0]++' > "$image"_TEMP_bottom
cat "$image"_TEMP_bottom | while read line;do cat RANGE | grep -w $line; done > "$image"_TEMP2_bottom
grep "\S" "$image"_TEMP2_bottom |  awk '!seen[$0]++' > "$image"_TEMP3_bottom
#cat RANGE_x | awk -v x="$StartX_bottom" -v xe="$EndX_bottom" '{if ($1>x && $1<xe) print}' > "$image"_TEMP_X_bottom

Acount_orange=0
Tcount_orange=0
Gcount_orange=0
Ccount_orange=0

Acount_orange_bottom=0
Tcount_orange_bottom=0
Gcount_orange_bottom=0
Ccount_orange_bottom=0


#HAPLOID REGION
if [ $Limit == 'NA' ] && [ $StartX != 'NA' ]
then
  Lines=$(cat "$image"_TEMP3)
  for line in $Lines
  do
    Yaxis=$(echo $line | awk '{print "y=""\""$1"\""}' | awk '{ gsub (" ", "", $0); print}')
    A=$(sed 's/<text/\n/g' $image | grep -w $Yaxis| awk '{gsub("x=\"","",$1);print}' | awk '{gsub("\"","",$1);print}' | awk -v x="$StartX" -v xe="$EndX" '{if ($1>x && $1<xe) print}' | grep -o '>A<' | wc -l)
    T=$(sed 's/<text/\n/g' $image | grep -w $Yaxis| awk '{gsub("x=\"","",$1);print}' | awk '{gsub("\"","",$1);print}' | awk -v x="$StartX" -v xe="$EndX" '{if ($1>x && $1<xe) print}' | grep -o '>T<' | wc -l)
    G=$(sed 's/<text/\n/g' $image | grep -w $Yaxis| awk '{gsub("x=\"","",$1);print}' | awk '{gsub("\"","",$1);print}' | awk -v x="$StartX" -v xe="$EndX" '{if ($1>x && $1<xe) print}' | grep -o '>G<' | wc -l)
    C=$(sed 's/<text/\n/g' $image | grep -w $Yaxis| awk '{gsub("x=\"","",$1);print}' | awk '{gsub("\"","",$1);print}' | awk -v x="$StartX" -v xe="$EndX" '{if ($1>x && $1<xe) print}' | grep -o '>C<' | wc -l)
    Acount_orange=$((Acount_orange+A))
    Tcount_orange=$((Tcount_orange+T))
    Gcount_orange=$((Gcount_orange+G))
    Ccount_orange=$((Ccount_orange+C)) 
    sum_tot=$((C+A+T+G))
    echo $sum_tot
    Total_orange=$((Acount_orange+Tcount_orange+Gcount_orange+Ccount_orange))
    Total_orange_bottom="NA"
  done >> "$image"_Sum_Line.txt
  Total_orange=$((Acount_orange+Tcount_orange+Gcount_orange+Ccount_orange))
  Total_orange_bottom="NA" 
  median=$(sort -n "$image"_Sum_Line.txt | awk ' { a[i++]=$1; } END { print a[int(i/2)]; }')
  median_bottom="NA"
  Total_read=$(wc -l "$image"_Sum_Line.txt |awk '{print $1}')
  Count_Orange_Upper_reads_with_interruptions=$(cat "$image"_Sum_Line.txt | awk '{if ($1 !=0) print $0}' | wc -l | awk '{print $1}')
  Total_read_bottom="NA"
  Total_orange_bottom="NA"
  Acount_orange_bottom="NA"
  Tcount_orange_bottom="NA"
  Gcount_orange_bottom="NA"
  Ccount_orange_bottom="NA"
  median_bottom="NA"
  Count_Orange_Bottom_reads_with_interruptions="NA"
fi

sum_tot=0

#DIPLOID REGION
if [ $Limit != 'NA' ] && [ $StartX != 'NA' ]
then
  Lines=$(cat "$image"_TEMP3)
  for line in $Lines
  do
    Yaxis=$(echo $line | awk '{print "y=""\""$1"\""}' | awk '{ gsub (" ", "", $0); print}')
    A=$(sed 's/<text/\n/g' $image | grep -w $Yaxis| awk '{gsub("x=\"","",$1);print}' | awk '{gsub("\"","",$1);print}' | awk -v x="$StartX" -v xe="$EndX" '{if ($1>x && $1<xe) print}' | grep -o '>A<' | wc -l)
    T=$(sed 's/<text/\n/g' $image | grep -w $Yaxis| awk '{gsub("x=\"","",$1);print}' | awk '{gsub("\"","",$1);print}' | awk -v x="$StartX" -v xe="$EndX" '{if ($1>x && $1<xe) print}' | grep -o '>T<' | wc -l)
    G=$(sed 's/<text/\n/g' $image | grep -w $Yaxis| awk '{gsub("x=\"","",$1);print}' | awk '{gsub("\"","",$1);print}' | awk -v x="$StartX" -v xe="$EndX" '{if ($1>x && $1<xe) print}' | grep -o '>G<' | wc -l)
    C=$(sed 's/<text/\n/g' $image | grep -w $Yaxis| awk '{gsub("x=\"","",$1);print}' | awk '{gsub("\"","",$1);print}' | awk -v x="$StartX" -v xe="$EndX" '{if ($1>x && $1<xe) print}' | grep -o '>C<' | wc -l)
    Acount_orange=$((Acount_orange+A))
    Tcount_orange=$((Tcount_orange+T))
    Gcount_orange=$((Gcount_orange+G))
    Ccount_orange=$((Ccount_orange+C))
    sum_tot=$((C+A+T+G))
    echo $sum_tot
  done >> "$image"_Sum_Line.txt
  Total_orange=$((Acount_orange+Tcount_orange+Gcount_orange+Ccount_orange))
  median=$(sort -n "$image"_Sum_Line.txt | awk ' { a[i++]=$1; } END { print a[int(i/2)]; }')
  Total_read=$(wc -l "$image"_Sum_Line.txt |awk '{print $1}') 
  Count_Orange_Upper_reads_with_interruptions=$(cat "$image"_Sum_Line.txt | awk '{if ($1 !=0) print $0}' | wc -l | awk '{print $1}')
fi
sum_tot_bottom=0

if [ $Limit != 'NA' ] && [ $StartX != 'NA' ]
then  
  Lines2=$(cat "$image"_TEMP3_bottom)
  for line2 in $Lines2
  do
    Yaxis_bottom=$(echo $line2 | awk '{print "y=""\""$1"\""}' | awk '{ gsub (" ", "", $0); print}')    
    A_bottom=$(sed 's/<text/\n/g' $image | grep -w $Yaxis_bottom| awk '{gsub("x=\"","",$1);print}' | awk '{gsub("\"","",$1);print}' | awk -v x="$StartX_bottom" -v xe="$EndX_bottom" '{if ($1>x && $1<xe) print}' | grep -o '>A<' | wc -l)
    T_bottom=$(sed 's/<text/\n/g' $image | grep -w $Yaxis_bottom| awk '{gsub("x=\"","",$1);print}' | awk '{gsub("\"","",$1);print}' | awk -v x="$StartX_bottom" -v xe="$EndX_bottom" '{if ($1>x && $1<xe) print}' | grep -o '>T<' | wc -l)
    G_bottom=$(sed 's/<text/\n/g' $image | grep -w $Yaxis_bottom| awk '{gsub("x=\"","",$1);print}' | awk '{gsub("\"","",$1);print}' | awk -v x="$StartX_bottom" -v xe="$EndX_bottom" '{if ($1>x && $1<xe) print}' | grep -o '>G<' | wc -l)
    C_bottom=$(sed 's/<text/\n/g' $image | grep -w $Yaxis_bottom| awk '{gsub("x=\"","",$1);print}' | awk '{gsub("\"","",$1);print}' | awk -v x="$StartX_bottom" -v xe="$EndX_bottom" '{if ($1>x && $1<xe) print}' | grep -o '>C<' | wc -l)
    Acount_orange_bottom=$((Acount_orange_bottom+A_bottom))
    Tcount_orange_bottom=$((Tcount_orange_bottom+T_bottom))
    Gcount_orange_bottom=$((Gcount_orange_bottom+G_bottom))
    Ccount_orange_bottom=$((Ccount_orange_bottom+C_bottom))      
    sum_tot_bottom=$((C_bottom+A_bottom+T_bottom+G_bottom))
    echo $sum_tot_bottom
  done >> "$image"_Sum_Line_bottom.txt
  median_bottom=$(sort -n "$image"_Sum_Line_bottom.txt | awk ' { a[i++]=$1; } END { print a[int(i/2)]; }')
  Total_read_bottom=$(wc -l "$image"_Sum_Line_bottom.txt |awk '{print $1}') 
  Total_orange_bottom=$((Acount_orange_bottom+Tcount_orange_bottom+Gcount_orange_bottom+Ccount_orange_bottom))
  Count_Orange_Bottom_reads_with_interruptions=$(cat "$image"_Sum_Line_bottom.txt | awk '{if ($1 !=0) print $0}' | wc -l | awk '{print $1}')
fi

if [ $Limit == 'NA' ] && [ $StartX == 'NA' ]
then
  Acount_orange="NA"
  Tcount_orange="NA"
  Gcount_orange="NA"
  Ccount_orange="NA"
  Total_orange="NA"
  median="NA"
  Total_read="NA"
  Total_orange_bottom="NA"
  Acount_orange_bottom="NA"
  Tcount_orange_bottom="NA"
  Gcount_orange_bottom="NA"
  Ccount_orange_bottom="NA"
  median_bottom="NA"
  Total_read_bottom="NA"
fi


rm "$image"_TEMP
rm "$image"_TEMP2
rm "$image"_TEMP3
rm "$image"_TEMP_bottom
rm "$image"_TEMP2_bottom
rm "$image"_TEMP3_bottom
rm RANGE_x
rm RANGE
rm "$image"_Sum_Line_bottom.txt
rm "$image"_Sum_Line.txt

printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' A_upper_Orange T_upper_Orange G_upper_Orange C_upper_Orange Total_upper_Orange Median_orange_upper Total_read_upper Count_Orange_Upper_reads_with_interruptions A_bottom_Orange T_bottom_Orange G_bottom_Orange C_bottom_Orange Total_bottom_Orange Median_orange_bottom Total_read_bottom Count_Orange_Bottom_reads_with_interruptions $Acount_orange $Tcount_orange $Gcount_orange $Ccount_orange $Total_orange $median $Total_read $Count_Orange_Upper_reads_with_interruptions $Acount_orange_bottom $Tcount_orange_bottom $Gcount_orange_bottom $Ccount_orange_bottom $Total_orange_bottom $median_bottom $Total_read_bottom $Count_Orange_Bottom_reads_with_interruptions
