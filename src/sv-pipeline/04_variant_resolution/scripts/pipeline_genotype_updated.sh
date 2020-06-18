set -e

wrkdir=/PHShome/hb875/SFARI/DeepSeq/HarrisonCompare/Batch_genotype/clean_genotype

cd $wrkdir

##Use SFARI Phase1 to prototype genotyper###


##Get breakpoint bed for SFARI variants##
svtk vcf2bed /data/talkowski/Samples/SFARI/deep_sv/asc_540/sv-pipeline/04_variant_resolution/ssc_20170911/batch_integration/SSC.vcf.gz SSC.int.bed -i SOURCES

###Generated a training file (1kg.train.loci.bed) of 64 multiallelic sites from 1kg Handsaker Study can be found in Rd test github repo### 



#########################################
#########################################
#########################################
                   RD
#########################################
#########################################
#########################################



##Generate Cutoffs using inital cutoffs of 0.25 (seed.txt) to classify variants classes##
covfile=/data/talkowski/Samples/SFARI/deep_sv/asc_540/bincov/matrices/ASC540.all.binCov.raw.bed.gz
medianfile=/data/talkowski/Samples/SFARI/deep_sv/asc_540/bincov/matrices/ASC540.all.binCov.median famfile=/data/talkowski/Samples/SFARI/lists/SFARI.519.fam
blacklist=/data/talkowski/Samples/SFARI/lists/SFARI.pilot.IDs.txt
whitelist=/data/talkowski/Samples/SFARI/lists/SFARI.batch1.IDs.txt
##create bed from 1kg for training##
awk -v id=$(head -n 1 $whitelist) '{if ($5=id) print}' /data/talkowski/hb875/code/RdTest_git/1kg.train.loci.bed|tr ' ' '\t'>$wrkdir/1kg.train.bed
bed="1kg.train.bed"
outputname="train"

Rscript /PHShome/hb875/talkowski/hb875/code/RdTest_git/RdTest.R -b $bed -c $covfile -m $medianfile -f $famfile  -n $outputname -l $blacklist -i 100000 -g TRUE -r /data/talkowski/hb875/code/RdTest_git/seed_cutoff.txt -y /data/talkowski/hb875/code/RdTest_git/bin_exclude.bed.gz

##Create cutoff file for depth##
Rscript /PHShome/hb875/talkowski/hb875/code/RdTest_git/generate_cutoff.R train.median_geno 4 1kgcutoff.txt

###Modify  copystate 1 & 3 1kgcutoff to account for sep cutoff discovered by random forest###
##Needs to come from a randomforest cut-off file ##
##Use CMC as example but manually modify file to actually match SFARI cutoffs##
cutoff_file=/data/talkowski/xuefang/data/CommonMind/sv-pipeline/03_variant_filtering_bak/cutoffs/CMC.cutoffs
sep=$( awk -F'\t' '{if ($1=="PESR" && $6==1000 && $5=="RD_Median_Separation") print $2}' $cutoff_file)
##get the highest depth dep##
sep_gt5kb=$( awk -F'\t' '{if ($1=="Depth" && $5=="RD_Median_Separation") print $2}' $cutoff_file|sort -nr|head -n 1)

##update cutoff file with new seps##
cat 1kgcutoff.txt| awk -v var=$sep '{if ($1=="1" && $4>1-var) $4=1-var; else if ($1=="2" && $4<1+var) $4=1+var; print}'|tr ' ' '\t'>pesrsepcutoff.txt

cat 1kgcutoff.txt| awk -v var=$sep_gt5kb '{if ($1=="1" && $4>1-var) $4=1-var; else if ($1=="2" && $4<1+var) $4=1+var; print}'|tr ' ' '\t'>depthsepcutoff.txt


##run across all Genotypes##
mkdir rdgeno_bed
cd rdgeno_bed

##use depth sep cutoff file for gt5kb so need to split separately###
less $wrkdir/SSC.int.bed|awk '{if (($5=="DEL" || $5=="DUP") && $3-$2>=5000) print $1,$2,$3,$4,$6,$5}'|fgrep -wf /data/talkowski/Samples/SFARI/lists/SFARI.batch1.IDs.txt|tr ' ' '\t'|split -l 100 -a 4 - bedgt5kb

##use pesrgt1kb sep cutoff file for lt1kb so need to split separately###
less $wrkdir/SSC.int.bed|awk '{if (($5=="DEL" || $5=="DUP") && $3-$2<5000) print $1,$2,$3,$4,$6,$5}'|fgrep -wf /data/talkowski/Samples/SFARI/lists/SFARI.batch1.IDs.txt|tr ' ' '\t'|split -l 100 -a 4 - bedpesr


##files for rdtest##
covfile=/data/talkowski/Samples/SFARI/deep_sv/asc_540/bincov/matrices/ASC540.all.binCov.raw.bed.gz
medianfile=/data/talkowski/Samples/SFARI/deep_sv/asc_540/bincov/matrices/ASC540.all.binCov.median famfile=/data/talkowski/Samples/SFARI/lists/SFARI.519.fam
blacklist=/data/talkowski/Samples/SFARI/lists/SFARI.pilot.IDs.txt

for file in bed*
do
bed=$file
outputname=$(echo $file|awk '{print $1".out"}')

if [[ $bed =~ "pesr" ]]; 
then
bsub -q short -o $bed.log "Rscript /PHShome/hb875/talkowski/hb875/code/RdTest_git/RdTest.R -b $bed -c $covfile -m $medianfile  -f $famfile  -n $outputname -l $blacklist -i 1000000 -g TRUE -r $wrkdir/pesrsepcutoff.txt -y /data/talkowski/hb875/code/RdTest_git/bin_exclude.bed.gz -v TRUE"
else
bsub -q short -o $bed.log "Rscript /PHShome/hb875/talkowski/hb875/code/RdTest_git/RdTest.R -b $bed -c $covfile -m $medianfile -f $famfile  -n $outputname -l $blacklist -i 1000000 -g TRUE -r $wrkdir/depthsepcutoff.txt -y /data/talkowski/hb875/code/RdTest_git/bin_exclude.bed.gz -v TRUE"
fi
done

##combine genotypes together###
cat *.geno|awk '!_[$0]++' >$wrkdir/rd.geno.all

##Per line sample##
rm rd.geno.cnv.bed
for column in {5..1920}
do
echo $column
id=$(head -n 1 rd.geno.all|awk -v var=$column '{print $var}')
awk -v var=$column -v var1=$id '{if ($var!=2) print $4 "\t" var1 "\t" $var}' rd.geno.all|fgrep -v cnvID >> rd.geno.cnv.bed

done

gzip rd.geno.cnv.bed


##QC##
##overall compare##
awk '{print $2}' pe.geno.freq.txt|cut -d"_" -f1-5|sort -u|fgrep -wf - SSC.int.bed|awk '{if ($NF!="depth") print $4"," $6}'|awk -F',' '{ for (i = 2; i <= NF; ++i)  print $1"_"$i  "\t" $i }'|fgrep -wf /data/talkowski/Samples/SFARI/lists/SFARI.batch1.IDs.txt|awk '{print $1}'>pe.int.ids.txt

zcat $wrkdir/pe.geno.final.txt.gz|awk '{if ($NF>0) print $1"_"$2}'>pe.final.ids.txt

##shared##
fgrep -wf pe.int.ids.txt pe.final.ids.txt|wc -l
##unique geno##
fgrep -wvf pe.int.ids.txt pe.final.ids.txt|wc -l
##unique initial##
fgrep -wvf pe.final.ids.txt pe.int.ids.txt|wc -l
##greater 5kb###
##shared##
fgrep -wf pe.int.ids.txt pe.final.ids.txt|egrep "DEL|DUP"|cut -d"_" -f1-5|fgrep -wf <(awk '{if ($3-$2>=5000) print $4}' cnv.bed)|wc -l
##unique geno##
fgrep -wvf pe.int.ids.txt pe.final.ids.txt|egrep "DEL|DUP"|cut -d"_" -f1-5|fgrep -wf <(awk '{if ($3-$2>=5000) print $4}' cnv.bed)|wc -l
##unique initial##
fgrep -wvf pe.final.ids.txt pe.int.ids.txt|egrep "DEL|DUP"|cut -d"_" -f1-5|fgrep -wf <(awk '{if ($3-$2>=5000) print $4}' cnv.bed)|wc -l




#########################################
#########################################
#########################################
                   PE
#########################################
#########################################
#########################################

##PE##

##Remove pilot only breakpoints##
##remove breakpoint at the same location###

zcat /data/talkowski/Samples/SFARI/deep_sv/asc_540/sv-pipeline/04_variant_resolution/ssc_20170911/batch_integration/SSC.vcf.gz| sed 's/0\/1/0\/0/g'|fgrep -wvf <(fgrep -wvf /data/talkowski/Samples/SFARI/lists/SFARI.batch1.IDs.txt  $wrkdir/SSC.int.bed|awk '{print $4}')  >pe.all.vcf


mkdir pegeno
cd $wrkdir/pegeno


##split vcf for efficiency##

cd /PHShome/hb875/SFARI/DeepSeq/HarrisonCompare/Batch_genotype/pe_chunk
egrep "^#" ../pe.all.vcf>pe.header.vcf

egrep -v "^#" $wrkdir/pe.all.vcf|split -l 10000 - 

petable=/data/talkowski/Samples/SFARI/deep_sv/asc_540/sv-pipeline/02_evidence_assessment/02b_petest/pe_counts/cohort.txt.gz
medianfile=/data/talkowski/Samples/SFARI/deep_sv/asc_540/bincov/matrices/ASC540.all.binCov.median 
whitelist=/data/talkowski/Samples/SFARI/lists/SFARI.batch1.IDs.txt

for i in x*
do
cat pe.header.vcf $i>pe.$i.vcf

bsub -q normal -o $i.log "svtk count-pe -s $whitelist --medianfile $medianfile pe.$i.vcf $petable $i.pe.out"
 
done

##Use PE as example but manually modify file to actually match SFARI cutoffs##
cutoff_file=/data/talkowski/xuefang/data/CommonMind/sv-pipeline/03_variant_filtering_bak/cutoffs/CMC.cutoffs
##convert PE p-value cutoff to count with poisson test assuming 0 background##
##I just manually change to SFARI since this will be much easier in python##
pe_count=$( awk -F'\t' '{if ( $5=="PE_log_pval") print $2}' $cutoff_file)

##SFARI cutoff##
pe_count=7


cat *.out|awk -v var=$pe_count '{if ($3>=var) print}'|fgrep -v name|gzip>$wrkdir/pe.geno.all.txt.gz

cd $wrkdir


##Use RD to find PE cutoffs##

###Develop some RD filters to ensure accurate PE genotypes##
##filter out potential overlapping variants##

egrep "DEL|DUP" SSC.int.bed>cnv.bed

bedtools intersect -wa -wb -a cnv.bed -b cnv.bed|awk '{if ($4 !=$11) print $4"\n" $11}'|sort -u>cnv.exclude.all.bed

##size filter (>1kb)##
awk '{if ($3-$2>=1000) print $4}' cnv.bed>size.pass.txt


##remove depth only##
awk '{if ($NF=="depth") print $4}' SSC.int.bed>depthonly.fail.txt

##only allow up to three copy state (0,1,2 or 1,2) for training##

tail -n +2 rd.geno.all|cut -f4-|fgrep -w 1|fgrep -w 2|awk '{ for (i = 2; i <= NF; ++i)  if($i>3) print $1 }'|sort -u>copystatefail.txt

tail -n +2 rd.geno.all|cut -f4-|fgrep -w 2|fgrep -w 3|awk '{ for (i = 2; i <= NF; ++i)  if($i>4 || $i<2) print $1 }'|sort -u>>copystatefail.txt

tail -n +2 rd.geno.all|cut -f4-|fgrep -wvf copystatefail.txt|fgrep -w 1|fgrep -w 2|awk '{print $1}'>copystate.pass.txt
tail -n +2 rd.geno.all|cut -f4-|fgrep -wvf copystatefail.txt|fgrep -w 2|fgrep -w 3|awk '{print $1}'>>copystate.pass.txt

##check for mulitcopystate##
awk '!seen[$1"_" $NF]++' rd.geno.cnv.bed|awk '{print $1}'|sort|uniq -c>cnv.copystate.multi.count.txt
less cnv.copystate.multi.count.txt|awk '{if ($1>3) print $2}'|fgrep -wf - cnv.bed|awk '{if ($3-$2<5000 && $3-$2>1000) print $1,$2,$3,$4,$6,$5}'|tr ' ' '\t'|egrep -v "^X|^Y">multi.test.bed

##remove repetitive ##
##breakpoint##
cat /data/talkowski/rlc47/src/GRch37.segdups_gaps_abParts_heterochrom.bed /data/talkowski/rlc47/src/GRCh37.*.RMSK.merged.bed|coverageBed -a <(awk '{print $1,$2,$2+1,$4 "\n" $1,$3,$3+1,$4}' cnv.bed|tr ' ' '\t')  -b -|awk '{if ($NF>0) print $4}'|sort -u>repeat.breakpoint.fail.ids.txt

##depth##
cat /data/talkowski/rlc47/src/GRCh37.*.RMSK.merged.bed /data/talkowski/rlc47/src/GRch37.segdups_gaps_abParts_heterochrom.lumpy.exclude.bed|sort -k1,1 -k2,2n |coverageBed -a <(sort -k1,1 -k2,2n cnv.bed) -b -|awk '{if ($NF>0.3) print $4}' >repeat.depth.fail.ids.txt

##Final filter##
##include for testing##
##remove X && Y###
awk '{if ($1!~"X" && $1!~"Y")print $4}' cnv.bed|fgrep -wf copystate.pass.txt|fgrep -wf size.pass.txt|fgrep -wvf <(cat cnv.exclude.all.bed depthonly.fail.txt  repeat.breakpoint.fail.ids.txt) >pe.train.include.txt



##Join RD and PE genotypes##

join -j 1  -a 1 -e "2" -o 1.2 1.3 1.4 2.4 <(zcat pe.geno.all.txt.gz|fgrep -wf pe.train.include.txt|awk '{print $1"_"$2 "\t" $0}'|sort -k1,1 ) <(fgrep -wf pe.train.include.txt rd.geno.cnv.bed|awk '{print $1"_"$2 "\t" $0}'|sort -k1,1)|tr ' ' '\t' > PE.RD.merged.txt 

##Get cutoffs to filter out incorrectly label hom in R#
Rscript generate_cutoff_PE.R <(awk '{print $1"_"$2"\t" $3 "\t" $4}' PE.RD.merged.txt) PE.all.cutoff.txt

##through out any copy state  calls that have reads less than with p=0.05 away from copy state 1 or 3##
het_cutoff=$(awk '{if ($1==3 || $1==1) print (1.645*$3) + $2 }' PE.all.cutoff.txt.metric)

##Rerun with excluding 0 or 4 copy states that fail het_cutoff###
awk -v var=$het_cutoff '{if (!(($4=="0" || $4=="4") && $3<var)) print}' PE.RD.merged.txt > PE.RD.hetfilter.merged.txt

##median##
##get median from copy state 0 and 4###
median_hom=$(awk '{if ($NF==0 || $NF==4) print $3}'  PE.RD.hetfilter.merged.txt|Rscript -e 'd<-scan("stdin", quiet=TRUE)' \
          -e 'median(d)'|tr '\n' '\t'|awk '{print $NF}')
##get std from 1 && 3  for hom restriction###          
sd_het=$(awk '{if ($NF==1 || $NF==3) print $3}'  PE.RD.hetfilter.merged.txt|Rscript -e 'd<-scan("stdin", quiet=TRUE)' \
          -e 'mad(d)'|tr '\n' '\t'|awk '{print $NF*1.645}')

##Go 1 SD outside of homozygous median##          

##Genotype PE genotype (0-ref, then estimate copy state based on copy state that is 1 sd from sd_het  )##
cat $wrkdir/pegeno/*.out|fgrep -v name|awk -v var=$pe_count -v var1=$median_hom -v var2=$sd_het '{if ($3<var) print $1,$2,$3,0;else print $1,$2,$3,int(($3/var1)+var2/var1)+1}'|gzip>$wrkdir/pe.geno.final.txt.gz


##QC check ##

##NOT necessary in pipeline###
##HWE##

##genotype frequency##
zcat $wrkdir/pe.geno.final.txt.gz|egrep "mo|fa"|awk '{if ($NF>1) print $1"_"2 ;else print $1"_"$NF}'|sort|uniq -c>pe.geno.freq.txt

awk '{print $2}' pe.geno.freq.txt|cut -d"_" -f1-5|sort -u|awk '{print 0 "\t" $1"_"0 "\n" 0" \t" $1"_"1 "\n" 0"\t" $1"_"2}' >blank.geno.txt

cat pe.geno.freq.txt blank.geno.txt|awk '!seen[$2]++'|awk '{print $0 "\t" $2}'|cut -d"_" -f1-10|sort -k2,2|awk '{if(a[$3])a[$3]=a[$3]"\t"$1; else a[$3]=$1;}END{for (i in a)print i "\t" a[i];}'|awk '{if ($2<958 && $1!~"X" && $1!~"Y")  print}'>pe.hwe.txt

cat pe.geno.freq.txt blank.geno.txt|awk '!seen[$2]++'|awk '{print $0 "\t" $2}'|cut -d"_" -f1-10|sort -k2,2|awk '{if(a[$3])a[$3]=a[$3]"\t"$1; else a[$3]=$1;}END{for (i in a)print i "\t" a[i];}'|awk '{if ($2<958 && $1!~"X" && $1!~"Y")  print}'|fgrep -wf size.pass.txt>pe.hwe.gt1kb.txt

##overall check##
awk '{print $2}' pe.geno.freq.txt|cut -d"_" -f1-5|sort -u|fgrep -wf - SSC.int.bed


##RD vs PE##
join -j 1  -a 1 -e "2" -o 1.1 1.2 2.2 <(zcat pe.geno.final.txt.gz|fgrep -wf pe.train.include.txt|awk '{print $1"_"$2 "\t" $NF}'|sort -k1,1 ) <(fgrep -wf pe.train.include.txt rd.geno.cnv.bed|awk '{print $1"_"$2 "\t" $NF}'|sort -k1,1)|tr ' ' '\t' > PE.RD.merged.assess.txt 

join -j 1  -a 1 -e "2" -o 1.1 1.2 2.2 <(zcat pe.geno.final.txt.gz|fgrep -wf pe.train.include.txt|awk '{print $1"_"$2 "\t" $3}'|sort -k1,1 ) <(fgrep -wf pe.train.include.txt rd.geno.cnv.bed|awk '{print $1"_"$2 "\t" $NF}'|sort -k1,1)|tr ' ' '\t' > PE.RD.merged.assess.count.txt 

 

less PE.RD.merged.assess.txt|awk '{if ($2!="0") print $2"_"$3}'|sort|uniq -c|sort -nrk1,1

##overall compare##
awk '{print $2}' pe.geno.freq.txt|cut -d"_" -f1-5|sort -u|fgrep -wf - SSC.int.bed|awk '{if ($NF!="depth") print $4"," $6}'|awk -F',' '{ for (i = 2; i <= NF; ++i)  print $1"_"$i  "\t" $i }'|fgrep -wf /data/talkowski/Samples/SFARI/lists/SFARI.batch1.IDs.txt|awk '{print $1}'>pe.int.ids.txt

zcat $wrkdir/pe.geno.final.txt.gz|awk '{if ($NF>0) print $1"_"$2}'>pe.final.ids.txt

##shared##
fgrep -wf pe.int.ids.txt pe.final.ids.txt|wc -l
##unique geno##
fgrep -wvf pe.int.ids.txt pe.final.ids.txt|wc -l
##unique initial##
fgrep -wvf pe.final.ids.txt pe.int.ids.txt|wc -l
##greater 5kb###
##shared##
fgrep -wf pe.int.ids.txt pe.final.ids.txt|egrep "DEL|DUP"|cut -d"_" -f1-5|fgrep -wf <(awk '{if ($3-$2>=5000) print $4}' cnv.bed)|wc -l
##unique geno##
fgrep -wvf pe.int.ids.txt pe.final.ids.txt|egrep "DEL|DUP"|cut -d"_" -f1-5|fgrep -wf <(awk '{if ($3-$2>=5000) print $4}' cnv.bed)|wc -l
##unique initial##
fgrep -wvf pe.final.ids.txt pe.int.ids.txt|egrep "DEL|DUP"|cut -d"_" -f1-5|fgrep -wf <(awk '{if ($3-$2>=5000) print $4}' cnv.bed)|wc -l

##get PE pass##
pe_cutoff=3.04


fgrep -wvf pe.final.ids.txt pe.int.ids.txt|egrep "DEL|DUP"|awk -F"_" '{print $0 "\t" $NF}'|awk '{print $0 "\t" $1}'|cut -d"_" -f1-10|awk '{print $3 "\t" $2}'|awk '{a[$1]=a[$1]?a[$1]","$NF:$NF;}END{for (i in a)print i "\t" a[i];}'|sort -k1,1|join -1 1 -2 4 - <(awk '{if ($3-$2>=5000) print}' cnv.bed|sort -k4,4)|awk '{print $3,$4,$5,$1,$2,$6}'|tr ' ' '\t' >missing.petest.bed


fgrep -wvf pe.final.ids.txt pe.int.ids.txt|egrep "DEL|DUP"|awk -F"_" '{print $0 "\t" $NF}'|awk '{print $0 "\t" $1}'|cut -d"_" -f1-10|awk '{print $3 "\t" $2}'|awk '{a[$1]=a[$1]?a[$1]","$NF:$NF;}END{for (i in a)print i "\t" a[i];}'|sort -k1,1|join -1 1 -2 4 - <(awk '{if ($3-$2>=5000) print}' cnv.bed|sort -k4,4)|awk '{print $3,$4,$5,$1,$2,$6}'|tr ' ' '\t' >missing.petest.bed


fgrep -wvf pe.final.ids.txt pe.int.ids.txt|egrep "DEL|DUP"|awk -F"_" '{print $0 "\t" $NF}'|awk '{print $0 "\t" $1}'|cut -d"_" -f1-10|fgrep -wf <(awk '{if ($3-$2>=5000) print $4}' cnv.bed)|awk '{print $3 "_" $2}'|fgrep -wf - <(zcat $wrkdir/pe.geno.final.txt.gz|awk '{print $1"_"$2 "\t" $3}')>missing.petest.count.txt


cat pe.int.ids.txt|awk -F"_" '{print $0 "\t" $NF}'|awk '{print $0 "\t" $1}'|cut -d"." -f1-3|egrep "mo|fa"|awk '{print $3}'|sort -u>missing.parent.pe.txt

fgrep -wvf pe.final.ids.txt pe.int.ids.txt|awk -F"_" '{print $0 "\t" $NF}'|awk '{print $0 "\t" $1}'|cut -d"." -f1-3|egrep "s1|p1"|fgrep -wf missing.parent.pe.txt|wc -l


fgrep -wvf pe.final.ids.txt pe.int.ids.txt|awk -F"_" '{print $0 "\t" $NF}'|awk '{print $0 "\t" $1}'|cut -d"." -f1-3|egrep "s1|p1"|fgrep -wvf missing.parent.pe.txt|wc -l

cat pe.final.ids.txt |awk -F"_" '{print $0 "\t" $NF}'|awk '{print $0 "\t" $1}'|cut -d"." -f1-3|egrep "mo|fa"|awk '{print $3}'|sort -u>missing.parent.pe.txt

fgrep -wvf pe.int.ids.txt pe.final.ids.txt |awk -F"_" '{print $0 "\t" $NF}'|awk '{print $0 "\t" $1}'|cut -d"." -f1-3|egrep "s1|p1"|fgrep -wf missing.parent.pe.txt|wc -l

fgrep -wvf pe.int.ids.txt pe.final.ids.txt |awk -F"_" '{print $0 "\t" $NF}'|awk '{print $0 "\t" $1}'|cut -d"." -f1-3|egrep "s1|p1"|fgrep -wvf missing.parent.pe.txt|wc -l


fgrep -wvf pe.int.ids.txt pe.final.ids.txt |egrep "DEL|DUP"|awk -F"_" '{print $0 "\t" $NF}'|awk '{print $0 "\t" $1}'|cut -d"_" -f1-10|fgrep -wf <(awk '{if ($3-$2>=5000) print $4}' cnv.bed)|egrep "s1|p1"|cut -d"." -f1|fgrep -wf missing.parent.pe.txt|wc -l

fgrep -wvf pe.int.ids.txt pe.final.ids.txt |egrep "DEL|DUP"|awk -F"_" '{print $0 "\t" $NF}'|awk '{print $0 "\t" $1}'|cut -d"_" -f1-10|fgrep -wf <(awk '{if ($3-$2>=5000) print $4}' cnv.bed)|egrep "s1|p1"|cut -d"." -f1|fgrep -wvf missing.parent.pe.txt|wc -l


fgrep -wvf pe.final.ids.txt pe.int.ids.txt |egrep "DEL|DUP"|awk -F"_" '{print $0 "\t" $NF}'|awk '{print $0 "\t" $1}'|cut -d"_" -f1-10|fgrep -wf <(awk '{if ($3-$2>=5000) print $4}' cnv.bed)|egrep "s1|p1"|cut -d"." -f1|fgrep -wf missing.parent.pe.txt|wc -l

fgrep -wvf pe.final.ids.txt pe.int.ids.txt|egrep "DEL|DUP"|awk -F"_" '{print $0 "\t" $NF}'|awk '{print $0 "\t" $1}'|cut -d"_" -f1-10|fgrep -wf <(awk '{if ($3-$2>=5000) print $4}' cnv.bed)|egrep "s1|p1"|cut -d"." -f1|fgrep -wvf missing.parent.pe.txt|wc -l


#########################################
#########################################
#########################################
                   SR
#########################################
#########################################
#########################################

##created SR correct vcf with proper coordinates for SR test (SR.vcf)##  
##Should not be an issue in future cohorts since Matt has added fix to now include these##
###get SR for all individuals##
##Pull out info to scan the SR matrix##
mkdir sr_geno


cd $wrkdir/sr_geno


egrep "^#" $wrkdir/SR.vcf>sr.header.vcf

egrep -v "^#" $wrkdir/SR.vcf|split -l 1000 - 

srtable=/data/talkowski/Samples/SFARI/deep_sv/asc_540/pesr_testing/split_counts/cohort.txt.gz
medianfile=/data/talkowski/Samples/SFARI/deep_sv/asc_540/bincov/matrices/ASC540.all.binCov.median

for i in x*
do
cat sr.header.vcf $i>sr.$i.vcf

bsub -q normal -o $i.log "svtk count-sr -s $whitelist --medianfile $medianfile sr.$i.vcf $srtable $i.sr.out"
 
done

cd $wrkdir

##get SR lower bound cut off###
cutoff_file=/data/talkowski/xuefang/data/CommonMind/sv-pipeline/03_variant_filtering_bak/cutoffs/CMC.cutoffs

##needs to be converted to count from poisson like pe###
sr_count=$( awk -F'\t' '{if ( $5=="SR_sum_log_pval") print $2}' $cutoff_file|head -n 1)

##set SR for SFARI##
sr_count=5


cat $wrkdir/sr_geno/*.out|fgrep -v name|gzip>combined.SR.normalized.counts.txt.gz


zcat combined.SR.normalized.counts.txt.gz|awk '{print $1"@"$3 "\t" $2 "\t" $4}' |awk '{a[$1]+=$3;}END{for(i in a)print i"\t"a[i];}' |tr '@' '\t'|sort -k1,1|gzip>combined.SR.normalized.sum.txt.gz

##Require both sides to have at least half of sr_count for training purposes##

zcat combined.SR.normalized.counts.txt.gz|awk -v sr_count=$sr_count '{if ($NF>(sr_count/2)) print $1"_"$3}'|sort|uniq -c|awk '{if ($1==2) print $2}'>two.sided.pass.txt

##Pass SR.test##
##need to cutoff results with SR test##
sr_cutoff=2.17147240951626

awk '{print $1}' ../combined.loc.ids.txt|fgrep -whf -  /data/talkowski/Samples/SFARI/deep_sv/asc_540/sv-pipeline/02_evidence_assessment/02c_srtest/srtest/*.stats|fgrep sum|awk -v var=$sr_cutoff '{if ($4>=var) print $1}'|fgrep -wf - ../combined.loc.ids.txt|awk '{print $2}'>pass.srtest.txt


###Test against RD##

##Join RD and SR genotypes and filter same as PE##
cat pe.train.include.txt>sr.train.include.txt

join -j 1  -a 1 -e "2" -o 1.2 1.3 1.4 2.4 <(zcat combined.SR.normalized.sum.txt.gz|fgrep -wf sr.train.include.txt|awk '{print $1"_"$2 "\t" $0}'|fgrep -wf two.sided.pass.txt|sort -k1,1 ) <(fgrep -wf sr.train.include.txt rd.geno.cnv.bed|awk '{print $1"_"$2 "\t" $0}'|fgrep -wf two.sided.pass.txt|sort -k1,1)|tr ' ' '\t' > SR.RD.merged.txt 

##Get cutoffs to filter out incorrectly label hom in R#
Rscript generate_cutoff_PE.R <(awk '{print $1"_"$2"\t" $3 "\t" $4}' SR.RD.merged.txt) SR.all.cutoff.txt

##through out any copy state  calls that have reads less than with p=0.05 away from copy state 1 or 3##
het_cutoff=$(awk '{if ($1==3 || $1==1) print (1.645*$3) + $2 }' SR.all.cutoff.txt.metric)

##Rerun without excluding 0 or 4 copy states that fail het_cutoff###
awk -v var=$het_cutoff '{if (!(($4=="0" || $4=="4") && $3<var)) print}' SR.RD.merged.txt > SR.RD.hetfilter.merged.txt

##median##
##get median from copy state 0 and 4###
median_hom=$(awk '{if ($NF==0 || $NF==4) print $3}'  SR.RD.hetfilter.merged.txt|Rscript -e 'd<-scan("stdin", quiet=TRUE)' \
          -e 'median(d)'|tr '\n' '\t'|awk '{print $NF}')
##get std from 1 && 3  for hom restriction###          
sd_het=$(awk '{if ($NF==1 || $NF==3) print $3}'  SR.RD.hetfilter.merged.txt|Rscript -e 'd<-scan("stdin", quiet=TRUE)' \
          -e 'mad(d)'|tr '\n' '\t'|awk '{print $NF*1.645}')

##Genotype SR genotype (0-ref, then estimate copy state based on copy state that is 1 sd from sd_het  )##



zcat combined.SR.normalized.sum.txt.gz|fgrep -wf pass.srtest.txt|awk '{print $0 "\t" $1"_"$2}'|fgrep -wf two.sided.pass.txt|cut -f1-3|awk -v var=$sr_count -v var1=$median_hom -v var2=$sd_het '{if ($3<var) print $1,$2,$3,0;else print $1,$2,$3,int(($3/var1)+var2/var1)+1}'>$wrkdir/sr.geno.final.txt

zcat combined.SR.normalized.sum.txt.gz|fgrep -wf pass.srtest.txt|awk '{print $0 "\t" $1"_"$2}'|fgrep -wvf two.sided.pass.txt|cut -f1-3|awk '{print $1,$2,$3,0}'>>$wrkdir/sr.geno.final.txt


gzip sr.geno.final.txt



zcat combined.SR.normalized.sum.txt.gz|fgrep -wf pass.srtest.txt|awk '{print $0 "\t" $1"_"$2}'|cut -f1-3|awk -v var=$sr_count -v var1=$median_hom -v var2=$sd_het '{if ($3<var) print $1,$2,$3,0;else print $1,$2,$3,int(($3/var1)+var2/var1)+1}'|gzip>$wrkdir/sr.geno.final.oneside.txt.gz



##can be difficult for interact node so sr genotype run on big##
bsub -q big -o sr.geno.log "./sr.geno.sh"



##QC##
##SR HWE##

##genotype frequency##
zcat $wrkdir/sr.geno.final.txt.gz|egrep "mo|fa"|awk '{if ($NF>1) print $1"_"2 ;else print $1"_"$NF}'|sort|uniq -c>sr.geno.freq.txt

awk '{print $2}' sr.geno.freq.txt|cut -d"_" -f1-5|sort -u|awk '{print 0 "\t" $1"_"0 "\n" 0" \t" $1"_"1 "\n" 0"\t" $1"_"2}' >blank.geno.txt

cat sr.geno.freq.txt blank.geno.txt|awk '!seen[$2]++'|awk '{print $0 "\t" $2}'|cut -d"_" -f1-10|sort -k2,2|awk '{if(a[$3])a[$3]=a[$3]"\t"$1; else a[$3]=$1;}END{for (i in a)print i "\t" a[i];}'|awk '{if ($2<958  && NF==4 && $1!~"X" && $1!~"Y")  print}'>sr.hwe.txt

cat sr.geno.freq.txt blank.geno.txt|awk '!seen[$2]++'|awk '{print $0 "\t" $2}'|cut -d"_" -f1-10|sort -k2,2|awk '{if(a[$3])a[$3]=a[$3]"\t"$1; else a[$3]=$1;}END{for (i in a)print i "\t" a[i];}'|awk '{if ($2<958 && NF==4 && $1!~"X" && $1!~"Y")  print}'|fgrep -wf size.pass.txt>sr.hwe.gt1kb.txt

join -j 1 <(sort -k1,1 ../cnv.size.txt) <(sort -k1,1 SR.RD.hetfilter.merged.txt)|awk '{print $2 "\t" $4 "\t" $5}'>SR_RD.plot.txt



##overall compare##
awk '{print $2}' sr.geno.freq.txt|cut -d"_" -f1-5|sort -u|fgrep -wf - SSC.int.bed|awk '{if ($NF!="depth") print $4"," $6}'|awk -F',' '{ for (i = 2; i <= NF; ++i)  print $1"_"$i  "\t" $i }'|fgrep -wf /data/talkowski/Samples/SFARI/lists/SFARI.batch1.IDs.txt|awk '{print $1}'>sr.int.ids.txt

zcat $wrkdir/sr.geno.final.txt.gz|awk '{if ($NF>0) print $1"_"$2}'>sr.final.ids.txt
zcat $wrkdir/sr.geno.final.oneside.txt.gz|awk '{if ($NF>0) print $1"_"$2}'>sr.final.ids.txt


##shared##
fgrep -wf sr.int.ids.txt sr.final.ids.txt|wc -l
##unique geno##
fgrep -wvf sr.int.ids.txt sr.final.ids.txt|wc -l
##unique initial##
fgrep -wvf sr.final.ids.txt sr.int.ids.txt|wc -l
##greater 5kb###
##shared##
fgrep -wf sr.int.ids.txt sr.final.ids.txt|egrep "DEL|DUP"|cut -d"_" -f1-5|fgrep -wf <(awk '{if ($3-$2>=5000) print $4}' cnv.bed)|wc -l
##unique geno##
fgrep -wvf sr.int.ids.txt sr.final.ids.txt|egrep "DEL|DUP"|cut -d"_" -f1-5|fgrep -wf <(awk '{if ($3-$2>=5000) print $4}' cnv.bed)|wc -l
##unique initial##
fgrep -wvf sr.final.ids.txt sr.int.ids.txt|egrep "DEL|DUP"|cut -d"_" -f1-5|fgrep -wf <(awk '{if ($3-$2>=5000) print $4}' cnv.bed)|wc -l



fgrep -wvf sr.final.ids.txt sr.int.ids.txt|awk -F"_" '{print $0 "\t" $NF}'|awk '{print $0 "\t" $1}'|cut -d"_" -f1-10|awk '{print $3 "_" $2}'|fgrep -wf - <(zcat $wrkdir/sr.geno.final.txt.gz|awk '{print $1"_"$2 "\t" $3}')>missing.srtest.count.txt

fgrep -wvf sr.final.ids.txt sr.int.ids.txt|awk -F"_" '{print $0 "\t" $NF}'|awk '{print $0 "\t" $1}'|cut -d"_" -f1-10|awk '{print $3 "_" $2}'|fgrep -wf - <(zcat $wrkdir/sr.geno.final.txt.gz|awk '{print $1"_"$2 "\t" $3}')>denovo.srtest.count.txt

cat sr.int.ids.txt|awk -F"_" '{print $0 "\t" $NF}'|awk '{print $0 "\t" $1}'|cut -d"." -f1-3|egrep "mo|fa"|awk '{print $3}'|sort -u>missing.parent.sr.txt

fgrep -wvf sr.final.ids.txt sr.int.ids.txt  |awk -F"_" '{print $0 "\t" $NF}'|awk '{print $0 "\t" $1}'|cut -d"." -f1-3|awk '{print $0 "\t" $NF}'|cut -d"_" -f1-15|egrep "DEL|DUP"|fgrep -wf <(awk '{if ($3-$2>=5000) print $4}' cnv.bed)|egrep "s1|p1"|fgrep -wf missing.parent.sr.txt|wc -l

fgrep -wvf sr.final.ids.txt sr.int.ids.txt  |awk -F"_" '{print $0 "\t" $NF}'|awk '{print $0 "\t" $1}'|cut -d"." -f1-3|awk '{print $0 "\t" $NF}'|cut -d"_" -f1-15|egrep "DEL|DUP"|fgrep -wf <(awk '{if ($3-$2>=5000) print $4}' cnv.bed)|egrep "s1|p1"|fgrep -wvf missing.parent.sr.txt|wc -l

cat sr.final.ids.txt |awk -F"_" '{print $0 "\t" $NF}'|awk '{print $0 "\t" $1}'|cut -d"." -f1-3|egrep "mo|fa"|awk '{print $3}'|sort -u>missing.parent.sr.txt

fgrep -wvf sr.int.ids.txt sr.final.ids.txt |awk -F"_" '{print $0 "\t" $NF}'|awk '{print $0 "\t" $1}'|cut -d"." -f1-3|awk '{print $0 "\t" $NF}'|cut -d"_" -f1-15|egrep "DEL|DUP"|fgrep -wf <(awk '{if ($3-$2>=5000) print $4}' cnv.bed)|egrep "s1|p1"|fgrep -wf missing.parent.sr.txt|wc -l

fgrep -wvf sr.int.ids.txt sr.final.ids.txt |awk -F"_" '{print $0 "\t" $NF}'|awk '{print $0 "\t" $1}'|cut -d"." -f1-3|awk '{print $0 "\t" $NF}'|cut -d"_" -f1-15|egrep "DEL|DUP"|fgrep -wf <(awk '{if ($3-$2>=5000) print $4}' cnv.bed)|egrep "s1|p1"|fgrep -wvf missing.parent.sr.txt|wc -l


##Find overlapping variants across all three day sets##
#join -j 1 -a 2 -e "0" -o 2.1 1.2 2.2 <(zcat sr.geno.final.txt.gz|awk '{if ($NF!=0) print $1"_"$2 "\t" $NF}'|sort -k1,1) <(zcat pe.geno.final.txt.gz|fgrep -wf <(cat SR1.vcf|egrep -v "^#"|awk '{print $3}')|awk '{if ($NF!=0) print $1"_"$2 "\t" $NF}'|sort -k1,1)|gzip>sr_pe.geno.combine.gz

#join -j 1 -a 2 -e "0" -o 2.1 1.2 2.2 <(zcat sr.geno.final.txt.gz|awk '{if ($NF!=0) print $1"_"$2 "\t" $NF}'|sort -k1,1) <(zcat pe.geno.final.txt.gz|fgrep -wf <(cat SR1.vcf|egrep -v "^#"|awk '{print $3}')|awk '{if ($NF!=0) print $1"_"$2 "\t" $NF}'|sort -k1,1)|gzip>sr_pe.geno.combine.gz


##join -j 1 -a 1 -e "0" -o 1.1 1.2 2.2 <(zcat sr.geno.final.txt.gz|fgrep -wf <(cat SR1.vcf|egrep -v "^#"|awk '{print $3}')|awk '{if ($NF!=0) print $1"_"$2 "\t" $NF}'|sort -k1,1) <(zcat pe.geno.final.txt.gz|awk '{if ($NF!=0) print $1"_"$2 "\t" $NF}'|sort -k1,1)|gzip>sr_pe.geno.combine.gz


cut -d"_" -f1-5 pe.geno.freq.txt|awk '{print $2}'|sort -u|fgrep -wf - <(cut -d"_" -f1-5 sr.geno.freq.txt|awk '{print $2}'|sort -u)>overlap.sr.pe.txt

##combine.sr.pe.sh##

join -j 1 -a 1 -e "0" -o 1.1 1.2 2.2 <(zcat sr.geno.final.txt.gz|awk '{ print $1"_"$2 "\t" $NF}'|sort -k1,1) <(zcat pe.geno.final.txt.gz|fgrep -wf overlap.sr.pe.txt|awk '{ print $1"_"$2 "\t" $NF}'|sort -k1,1)|gzip>sr_pe.geno.combine.gz

bsub -q big -o srpe.combine.log "./combine.sr.pe.sh"

##remove X && Y because we are not accounting for sex yet##
join -j 1 -a 1 -e "2" -o 1.1 1.2 1.3 2.2 <(zcat sr_pe.geno.combine.gz) <(cat rd.geno.cnv.bed|fgrep -wf overlap.sr.pe.txt|awk '{ print $1"_"$2 "\t" $NF}'|sort -k1,1)|awk '{if ($1!~"X" && $1!~"Y") print}'|gzip>sr_pe_rd.geno.combine.gz

zcat sr_pe_rd.geno.combine.gz|awk '{if (!($2==0 && $3==0 && $4==2)) print}'|gzip>sr_pe_rd.geno.noref.combine.gz

zcat sr_pe_rd.geno.noref.combine.gz|awk '{if ($2>3) print $1,3,$3,$4;else print $0}'|awk '{if ($3>3) print $1,$2,3,$4;else print $0}'|awk '{if ($4>5) print $1,$2,$3,5;else print $0 }'|gzip>sr_pe_rd.geno.noref.nomulti.combine.gz

##Rd multi##
zcat  sr_pe_rd.geno.noref.nomulti.combine.gz|awk '{print $0 "\t" $1}'|cut -d"_" -f1-10|awk '!seen[$4"_"$5]++'|awk '{print $NF}'|sort|uniq -c|awk '{if ($1>3) print $2}'>multi.ids.txt

zcat  sr_pe_rd.geno.noref.nomulti.combine.gz|awk '{print $2"_"$3"_"$4}'|sort|uniq -c|sort -nk1,1|awk '{print $1"_"$2}'

##gt 5kb##
zcat  sr_pe_rd.geno.noref.nomulti.combine.gz|awk '{print $0 "\t" $1}'|cut -d"_" -f1-10|fgrep -wf <(awk '{if ($3-$2>=5000) print $4}' cnv.bed)|fgrep -wvf multi.ids.txt|awk '{print $2"_"$3"_"$4}'|sort|uniq -c|sort -nk1,1|awk '{print $1"_"$2}'


zcat  sr_pe_rd.geno.noref.nomulti.combine.gz|awk '{print $0 "\t" $1}'|cut -d"_" -f1-10|fgrep -wf <(awk '{if ($3-$2>=5000) print $4}' cnv.bed)|fgrep -wvf multi.ids.txt|fgrep -wf pass.srtest.txt|awk '{if($2==0 && $3==0) print $0 "\t" $1}'|awk -F"_" '{print $0 "\t" $NF}'|awk '{print $5 "\t" $NF}'|awk '{a[$1]=a[$1]?a[$1]","$2:$2;}END{for (i in a)print i "\t" a[i];}'|sort -k1,1|join -1 1 -2 4 - <(sort -k4,4 cnv.bed)|awk '{print $3,$4,$5,$1,$2,$6}'|tr ' ' '\t'>test.bed



zcat sr_pe.geno.combine.gz|awk -F'_' '{print $4 "\t" $NF}'|awk '{print $1"_" $3 "_" $4 }'|sort|uniq -c|sort -nk1,1

zcat sr_pe.geno.combine.gz|awk '{if ($2=="0" && $3>0) print }'|sort -k1,1|join -a 1 -e '2' -j 1 -o 1.1 1.3 2.2 - <(awk '{print $1"_"$2 "\t" $3}' rd.geno.cnv.bed|sort -k1,1)>peonly.rd.txt

cat peonly.rd.txt|awk -F'_' '{print $4 "\t" $NF}'|awk '{print $1"_" $4 "_" $5 }'|sort|uniq -c|sort -nk1,1


zcat sr_pe.geno.combine.gz|awk '{if ($2=="1" && $3=="1") print }'|cut -d"_" -f1-5|sort -u|fgrep -wvf - <(cat peonly.rd.txt|awk '{if ($2=="1" && $3=="1") print $0 "\t" $1}'||fgrep DEL|cut -d"_" -f1-10)|awk '{print $NF "\t" $1}'|awk -F"_" '{print $0 "\t" $NF}'|awk '{a[$1]=a[$1]?a[$1]","$NF:$NF;}END{for (i in a)print i "\t" a[i];}'|sort -k1,1|join -1 1 -2 4 - <(sort -k4,4 cnv.bed) |awk '{print $3,$4,$5,$1,$2,$6}'|tr ' ' '\t'|awk '{if ($3-$2<1000) print}'|sort -R|head -n 50 >test.bed


zcat sr_pe.geno.combine.gz|awk '{if ($2=="0" && $3=="1") print }'|cut -d"_" -f1-5|sort -u|fgrep -wf - <(cat peonly.rd.txt|awk '{if ($2=="1" && $3=="2") print $0 "\t" $1}'|fgrep DEL|cut -d"_" -f1-10)|awk '{print $NF "\t" $1}'|awk -F"_" '{print $0 "\t" $NF}'|awk '{a[$1]=a[$1]?a[$1]","$NF:$NF;}END{for (i in a)print i "\t" a[i];}'|sort -k1,1|join -1 1 -2 4 - <(sort -k4,4 cnv.bed) |awk '{print $3,$4,$5,$1,$2,$6}'|tr ' ' '\t'|awk '{if ($3-$2<1000) print}'|egrep -v "^X|^Y" >nosr_rd_pe.bed



##include repeat regions in test##
awk '{if ($1!~"X" && $1!~"Y")print $4}' cnv.bed|fgrep -wf <(awk '{if ($3-$2>500) print $4}' cnv.bed)|fgrep -wvf <(cat cnv.exclude.all.bed depthonly.fail.txt) >qc.train.include.txt

join -j 1 -a 2 -e "0" -o 2.1 1.2 2.2 <(zcat sr.geno.final.txt.gz|awk '{if ($NF!=0) print $1"_"$2 "\t" $NF}'|sort -k1,1) <(cat rd.geno.cnv.bed|fgrep -wf qc.train.include.txt|fgrep -wf <(cat SR.vcf|egrep -v "^#"|awk '{print $3}')|awk '{if ($NF!=2) print $1"_"$2 "\t" $NF}'|sort -k1,1)|gzip>sr_rd.geno.combine.gz

##test how many sr missing rd support##
join -j 1 -a 1 -e "2" -o 1.1 1.2 2.2 <(zcat sr.geno.final.txt.gz|fgrep -wf qc.train.include.txt|awk '{if ($NF!=0) print $1"_"$2 "\t" $NF}'|sort -k1,1) <(cat rd.geno.cnv.bed|awk '{if ($NF!=2) print $1"_"$2 "\t" $NF}'|sort -k1,1)|gzip>rd_sr.geno.combine.gz

##plot by size##
see R code from ##~/Documents/Code/batch_genotype_updated_3_20_18.sh##

zcat sr_rd.geno.combine.gz|awk -F'_' '{print $4 "\t" $NF}'|awk '{print $1"_" $3 "_" $4 }'|sort|uniq -c|sort -nk1,1


join -j 1 -a 2 -e "0" -o 2.1 1.2 2.2 <(zcat pe.geno.final.txt.gz|awk '{if ($NF!=0) print $1"_"$2 "\t" $NF}'|sort -k1,1) <(cat rd.geno.cnv.bed|fgrep -wf qc.train.include.txt|awk '{if ($NF!=2) print $1"_"$2 "\t" $NF}'|sort -k1,1)|gzip>pe_rd.geno.combine.gz

zcat pe_rd.geno.combine.gz|awk -F'_' '{print $4 "\t" $NF}'|awk '{print $1"_" $3 "_" $4 }'|sort|uniq -c|sort -nk1,1

##PE_RD##
join -j 1 -a 1 -e "2" -o 1.1 1.2 2.2 <(zcat pe.geno.final.txt.gz|fgrep -wf qc.train.include.txt|awk '{if ($NF!=0) print $1"_"$2 "\t" $NF}'|sort -k1,1) <(cat rd.geno.cnv.bed|awk '{if ($NF!=2) print $1"_"$2 "\t" $NF}'|sort -k1,1)|gzip>rd_pe.geno.combine.gz

##SR before vs after##
egrep -v "^#" $wrkdir/SR.vcf|awk '{print $3}'|sort -u>variants.withsr.test

zcat combined.SR.normalized.sum.txt.gz|awk '{print $1}'|awk '!seen[$0]++'>somecount.txt

fgrep -wf pass.srtest.txt cnv.bed|fgrep -wf somecount.txt|fgrep -wf variants.withsr.test|cut -f1-6|awk -F'[,\t]' '{for (i=6;i<=NF;i++) print $4"_"$i }'|fgrep -f /data/talkowski/Samples/SFARI/lists/SFARI.batch1.IDs.txt>SR.inital.txt

fgrep -wvf two.sided.pass.txt SR.inital.txt|fgrep -wf <(zcat combined.SR.normalized.sum.txt.gz|awk '{if ($3>4) print $1"_"$2}')|awk '{print $0 "\t" $1}'|cut -d"_" -f1-10|less

fgrep -wf SR.inital.txt two.sided.pass.txt|fgrep -wf - <(zcat sr.geno.final.txt.gz|awk '{if ($NF!=0) print $1"_"$2 "\t" $0 }') |fgrep -wvf size.pass.txt|awk '{print $4}'|sort|uniq -c|sort -nk1,1

##de novo check##
while read line
do

echo $line|awk '{print $6}' |tr ',' '\n'|egrep "mo|fa"|awk -F'.' '{print $1}'>parent.txt
echo $line|awk '{print $6}' |tr ',' '\n'|egrep "s1|p1">child.txt

denovo_rate=$(paste <(fgrep -vf parent.txt child.txt|wc -l ) <(cat child.txt|wc -l)|awk '{print $1/$2}' )

echo $line|awk -v var=$denovo_rate '{print $4 "\t" var  }' >>denovo_rate.txt

rm parent.txt
rm child.txt

done< <(egrep "s1|p1" cnv.bed)


fgrep -v "END=;" SR.vcf|egrep -v "#"|awk -F'[\t;]' '{print $2"_"$8 "\t" $0}'|sed 's/_END=/\t/g'|awk '{if ($2-$1>=50) print}'|cut -f3-|cat <(egrep ^"#" SR.vcf) - >SR1.vcf

