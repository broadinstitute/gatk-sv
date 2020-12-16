BEGIN {FS=OFS="\t"}

# initialize a map of integer EV values onto string EVs
BEGIN {EVval[1]="RD"; EVval[2]="PE"; EVval[3]="RD,PE"; EVval[4]="SR";
	EVval[5]="RD,SR"; EVval[6]="PE,SR"; EVval[7]="RD,PE,SR";}

# read the names of the allosomes from an fai file
BEGIN {while ( getline < allosomeFile )
	{if ( $1~/X/ ) xChr = $1; else if ( $1~/Y/ ) yChr = $1}
	 if ( xChr=="" )
	 	{print "Can't determine the name of the X chromosome from " allsomeFile > "/dev/stderr"; exit 1}
	 if ( yChr=="" )
	 	{print "Can't determine the name of the Y chromosome from " allsomeFile > "/dev/stderr"; exit 1}}

# read the sexes of the samples from the ped file
BEGIN {while ( getline < pedFile ) sexForSampleName[$2] = $5}

# read the background file
BEGIN {while ( getline < bgdFile ) bgdEvent[$1]}

# read the both-sides file
BEGIN {while ( getline < bothFile ) bothEvent[$NF]}

# regexp for the #CHROM line with all the sample names
/^#C/  {# add filters at the end of the metadata
	print "##FILTER=<ID=HIGH_SR_BACKGROUND,Description=\"High number of SR splits in background samples indicating messy region\">";
	print "##FILTER=<ID=UNRESOLVED,Description=\"Variant is unresolved\">";
	print "##FILTER=<ID=BOTHSIDES_SUPPORT,Description=\"Variant has read-level support for both sides of breakpoint\">";
	# print the list of samples to includelist.txt
	for ( i=10; i<=NF; i++ )
		{print $i > "includelist.txt";
		 if ( !($i in sexForSampleName) )
		 	{print "Can't find sex for " $i " in " pedFile > "/dev/stderr"; exit 1}
		 sampleSex[i] = sexForSampleName[$i];}}

# fix the type of the EV field in the FORMAT line for EV
/^##FORMAT=<ID=EV,/ {$0="##FORMAT=<ID=EV,Number=.,Type=String,Description=\"Classes of evidence supporting final genotype\">"}

# regexp for all the non-header lines
!/^#/  {# replace integer EVs with strings
	nFmts=split($9,fmt,":");
	for ( fmtIdx=1; fmtIdx<=nFmts; fmtIdx++ )
		if ( fmt[fmtIdx]=="EV" ) break;
	if ( fmtIdx<=nFmts )   # if we found the EV field
		for ( sampIdx=10; sampIdx<=NF; ++sampIdx )   # for each sample
			# split out the genotype fields and replace the EV value
	 		{split($sampIdx,gt,":"); gt[fmtIdx]=EVval[gt[fmtIdx]];
			 # re-join the genotype fields for this sample
			 $sampIdx=gt[1];
			 for ( gtIdx=2; gtIdx<=nFmts; ++gtIdx )
				 $sampIdx=$sampIdx ":" gt[gtIdx]}

	# do some processing on the INFO fields
	nInfos=split($8,info,";");
	$8=""; sep=""; svType=""; end=0;
	for ( infoIdx=1; infoIdx<=nInfos; ++infoIdx )
		{infoFld=info[infoIdx];
		 # replace the alt allele with the SVTYPE info field, except for MEs
		 if ( infoFld~/^SVTYPE=/ )
		 	{svType = substr(infoFld,8);
			 if ( $5!~/:ME/ ) $5="<" svType ">"}
		 else if ( infoFld~/^END=/ )
		 	 end=substr(infoFld,5);
		 # set $6 to varGQ and remove varGQ field
		 else if ( infoFld~/^varGQ=/ ) {$6=substr(infoFld,7); continue}
		 # remove UNRESOLVED flags from INFO and move them to FILTERs
		 else if ( infoFld=="UNRESOLVED" ) {$7=$7 ";UNRESOLVED"; continue}
		 # remove MULTIALLELIC flags from INFOs
		 else if ( infoFld=="MULTIALLELIC" ) continue;
		 $8=$8 sep infoFld; sep=";"}

	# mark the background and both-sides events
	if ( $3 in bgdEvent ) $7=$7 ";HIGH_SR_BACKGROUND";
	if ( $3 in bothEvent ) {if ( $7=="PASS" ) $7="BOTHSIDES_SUPPORT"; else $7=$7 ";BOTHSIDES_SUPPORT"}

	# patch some genotypes on allosomes
	if ( $1==xChr || $1==yChr )
		{# non-male, non-female sex gets a GT of "./."
		 for ( sampIdx=10; sampIdx<=NF; ++sampIdx )
			{sampSex=sampleSex[sampIdx];
			 if ( sampSex!=1 && sampSex!=2 ) sub(/[^;]*;/,"./.",$i)}
		 # males are revised only for certain large DELS and DUPS that have specific median(RD_CN) values for each sex
		 if ( (svType=="DEL" || svType=="DUP") && end-$2+1 > 5000 )
			{for ( fmtIdx=1; fmtIdx<=nFmts; ++fmtIdx )
				if ( fmt[fmtIdx]=="RD_CN" ) break;
			 if ( fmtIdx<=nFmts ) # if we found the RD_CN field
			 	{male[0] = male[1] = male[2] = male[3] = 0;
				 female[0] = female[1] = female[2] = female[3] = 0;
			 	 for ( sampIdx=10; sampIdx<=NF; ++sampIdx )
					{split($sampIdx,gt,":");
					 rdCN=gt[fmtIdx];
					 if ( rdCN=="." ) continue;
					 if ( rdCN>2 ) rdCN = 3;
					 sampSex=sampleSex[sampIdx];
					 if ( sampSex==1 ) ++male[rdCN];
					 else if ( sampSex==2 ) ++female[rdCN];}
				 maleMedianCount=(male[0]+male[1]+male[2]+male[3])/2;
				 counts = 0;
				 for ( rdcnIdx=0; rdcnIdx<=3; ++rdcnIdx )
					{counts += male[rdcnIdx];
					 if ( counts==maleMedianCount ) {maleMedian = rdcnIdx + .5; break;}
					 else if ( counts>maleMedianCount ) {maleMedian = rdcnIdx; break;}}
				 counts = 0;
				 femaleMedianCount=(female[0]+female[1]+female[2]+female[3])/2;
				 for ( rdcnIdx=0; rdcnIdx<=3; ++rdcnIdx )
					{counts += female[rdcnIdx];
					 if ( counts==femaleMedianCount ) {femaleMedian = rdcnIdx + .5; break;}
					 else if ( counts>femaleMedianCount ) {femaleMedian = rdcnIdx; break;}}
				 if ( $1==xChr && maleMedian==1 && femaleMedian==2 ||
				      $1==yChr && maleMedian==1 && femaleMedian==0 )
				 	{print $3 > "sexchr.revise.txt";
					 for ( sampIdx = 10; sampIdx<=NF; ++sampIdx )
						{if ( sampleSex[sampIdx]!=1 ) continue;
						 split($sampIdx,gt,":");
						 rdCN=gt[fmtIdx];
					 	 if ( svType=="DEL" )
						 	{if ( rdCN>=1 ) gt[1]="0/0";
							 else if ( rdCN==0 ) gt[1] = "0/1"}
						 else # svType is "DUP"
						 	{if ( rdCN<=1 ) gt[1]="0/0";
							 else if ( rdCN==2 ) gt[1]="0/1";
							 else gt[1]="1/1"}
						 gt[fmtIdx]=rdCN + 1;
						 $sampIdx=""; sep="";
						 for ( gtIdx=1; gtIdx<=nFmts; ++gtIdx )
						 	{$sampIdx=$sampIdx sep gt[gtIdx]; sep=":";}}}}}

		 # females have all formatted fields erased for Y events
		 if ( $1==yChr )
		 	{emptyGT = "./.";
			 for ( fmtIdx=2; fmtIdx<=nFmts; ++fmtIdx )
				emptyGT=emptyGT ":.";
			 for ( sampIdx=10; sampIdx<=NF; ++sampIdx )
				if ( sampleSex[sampIdx]==2 ) $sampIdx=emptyGT}
		} #end of "if allosome"

	# ref allele is always "N"
	$4="N"}

{print $0}

