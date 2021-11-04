# WGD
**W**hole-**G**enome **D**osage: a suite of tools to evaluate dosage in whole-genome sequencing libraries

---  
## Table of Contents  
#### Examples  
- [Example WGD workflow](https://github.com/RCollins13/WGD#example-wgd-workflow)  
- [Pipeline runtimes](https://github.com/RCollins13/WGD#pipeline-runtimes)
- [CNV visualization and annotation](https://github.com/RCollins13/WGD#cnv-visualization-and-annotation)

#### Script documentation  
- [binCov.py](https://github.com/RCollins13/WGD#bincovpy)  
- [WG_binCov.sh](https://github.com/RCollins13/WGD#wg_bincovsh)  
- [makeMatrix.sh](https://github.com/RCollins13/WGD#makematrixsh)  
- [medianCoverage.R](https://github.com/RCollins13/WGD#mediancoverager)
- [compressCov.sh](https://github.com/RCollins13/WGD#compresscovsh)  

--- 

### Example WGD workflow  
#### Important Note Regarding Bedtools Compatibility  
As of ```bedtools 2.24.0```, some arguments and parameters changed for ```bedtools coverage```, an essential element of the WGD workflow. The WGD pipeline has been tested with ```bedtools 2.20.1``` and ```bedtools 2.25.1```, thus includes compatibility for both "new" and "old" parameter specifications, but one step (```binCov.py```) will require an argument to specify which version of ```bedtools``` you are using. See more details in the [usage note for ```binCov.py```](https://github.com/RCollins13/WGD#bincovpy) 
#### Prerequisites  
The WGD pipeline requires the following:  
- Coordinate-sorted, indexed bams for all samples
- List of contigs to evaluate
- Bed-file of N-masked regions of the reference genome. These are available for most reference genome assemblies from [UCSC](http://genome.ucsc.edu/ "UCSC Genome Browser")  
- 6F adjustment metadata matrices, which are available for human references hg19 or GRCh37 upon request. Matrices for other references are forthcoming.  

#### Step 1: Generate binned coverage per chromosome on all libraries  
Binned coverage is calculated by ```binCov.py``` on a per-chromosome basis. For human whole-genome sequencing analyses, nucleotide coverage at bin sizes of 100bp is recommended. Only primary contigs recommended; e.g. 1...Y. Parallelization of this process is encouraged.  

There are two approaches to parallelization, depending on your available computational resources. Examples are given below using LSF as a scheduler, but could be easily configured to your scheduler/environment.  

**Fully parallelized approach** (almost always preferable)
```
#!/usr/env/bin bash  
while read sample; do
  while read contig; do
    bsub "binCov.py ${sample}.bam ${contig} ${sample}.${contig}.rawCov.bed \
    -n ${sample}.${contig}.normCov.bed \
    -b 100 \
    -m nucleotide \
    -x /path/to/Nmask.bed"
  done < list_of_contigs.txt
done < list_of_samples.txt
```  

Alternatively, if available cores are limited or sample sizes are large, a wrapper script, ```WG_binCov.sh```, will calculate normalized coverage for a set of contigs in serial from a single bam.  

**Semi-parallelized approach** (preferred if available cores are dramatically fewer than [#contigs x #samples])
```
#!/usr/env/bin bash  
while read sample; do
  bsub "WG_binCov.sh ${sample}.bam ${sample} `pwd` \
  -L list_of_contigs.txt \
  -b 100 \
  -m nucleotide \
  -x /path/to/Nmask.bed"
done < list_of_samples.txt
```

#### Step 2: 6F coverage correction
Once you have 100bp binCov files per chromosome computed for all samples, the next step is to apply 6F correction. This is accomplished with the script ```multiCorrection.R``` in ```WGD/bin/```.  

To run 6F correction, the ```multiCorrection.R``` script requires a tab-delimited input file with the following three columns for each chromosome:  
1. full path to raw binCov file  
2. desired full path to output corrected binCov file  
3. full path to 6F metadata matrix for that chromosome  

The latter component, the 6F metadata matrix, is currently available upon request for hg19 and GRCh37. Additional references are forthcoming.  

An example head & tail of the expected input file for `multiCorrection.R` should be:  
```
/path/to/raw_binCov.chr1.bed.gz /path/to/output_corrected_binCov.chr1.bed.gz  /path/to/6F_correction_matrix.chr1.bed.gz
/path/to/raw_binCov.chr2.bed.gz /path/to/output_corrected_binCov.chr2.bed.gz  /path/to/6F_correction_matrix.chr2.bed.gz
/path/to/raw_binCov.chr3.bed.gz /path/to/output_corrected_binCov.chr3.bed.gz  /path/to/6F_correction_matrix.chr3.bed.gz
...
/path/to/raw_binCov.chr22.bed.gz  /path/to/output_corrected_binCov.chr22.bed.gz /path/to/6F_correction_matrix.chr22.bed.gz
/path/to/raw_binCov.chrX.bed.gz /path/to/output_corrected_binCov.chrX.bed.gz  /path/to/6F_correction_matrix.chrX.bed.gz
/path/to/raw_binCov.chrY.bed.gz /path/to/output_corrected_binCov.chrY.bed.gz  /path/to/6F_correction_matrix.chrY.bed.gz
```

Once you've compiled this input file per sample, you can call the 6F correction script as follows:  
```
$ WGD/bin/multiCorrection.R -z /path/to/6F_correction_input_file.txt
```

This job should require less than a couple hours to run on most standard machines.  


#### Step 3: Run WGD dosage scoring model  
Once 6F correction is complete for all chromosomes per sample, you next need to run the WGD dosage bias scoring model. This has three intermediate steps, detailed below.  

**Step 2.1:** Generate a WGD scoring input file per sample. This is done by combining all 6F-adjusted binCov files across all chromosomes, then filtering them down to a specific set of bins we know to be informative for dosage bias scoring. This is accomplished with a four-liner, as follows:
```
zcat /path/to/output_corrected_binCov.chr*.bed.gz | \
  bedtools intersect -f 1.0 -r -wa -a - \
  -b WGD/refs/WGD_scoring_mask.6F_adjusted.100bp.h37.bed > \
  /path/to/WGD_scoring_file_output.bed
```

**Step 2.2:** Combine WGD scoring input files across all samples into a WGD scoring matrix. This is done the same way you would create any binCov matrix (e.g. before running cn.MOPS or RdTest), using the ```makeMatrix.sh``` script in ```WGD/bin/```. Example usage:
```
#First, generate makeMatrix input file
while read ID; do
  echo -e "${ID}\t/path/to/${ID}.6F_adjusted.WGD_scoring_bins.bed.gz"
done < samples.list > \
/path/to/WGD_scoring.makeMatrix_input.txt

#Next, make the matrix
/WGD/bin/makeMatrix.sh -z \
    -o /path/to/output_WGD_scoring_masked.matrix.bed \
    -r /WGD/refs/WGD_scoring_mask.6F_adjusted.100bp.h37.bed \
    /path/to/WGD_scoring.makeMatrix_input.txt
```

**Step 2.3:** Run the WGD scoring model. Finally, you run the WGD scoring model on the matrix you created above with the following command:
```
/WGD/bin/scoreDosageBiases.R -z \
  -O /path/to/WGD_model_output/ \
  /path/to/output_WGD_scoring_masked.matrix.bed.gz \
  /WGD/refs/WGD_scoring_mask.6F_adjusted.100bp.h37.bed
```

The output from the model will be a list of dosage scores, and a pdf of two plots visualizing the dosage score distributions for all samples.

#### Step 4: Run ploidy estimation model
Similar to the WGD scoring model above, the ploidy estimation model has four sequential steps, detailed below.  

**Step 3.1:** Recompress binCov data per sample per chromosome to 1Mb bins. This is accomplished using the ```compressCov.sh``` script shipped with the WGD repo (under ```WGD/bin/```). Example usage below:  
```
#Do this per chromosome per sample
/WGD/bin/compressCov.sh -N -z -s \
    -o /path/to/6F_adjusted.1Mb_bins.binCov.bed.gz \
    /path/to/6F_adjusted.100bp_bins.binCov.bed.gz \
    10000
```

**Step 3.2:** Combine all recompressed binCov files into a single file per sample. This is done using the output from step 3.1, and can be accomplized using something simple like ```zcat``` or ```cat```.  

**Step 3.3:** Create matrix of 1Mb recompressed binCov values across all samples. This is done using the output from step 3.2, and can be done using the same instructions as in step 2.2.  

**Step 3.4:** Run the ploidy estimation model. This is performed on the 1Mb coverage matrix from step 3.3. with the ```estimatePloidy.R``` script. Example usage below:  
```
/WGD/bin/estimatePloidy.R -z \
    -O /path/to/ploidy_estimate_output_directory/ \
    /path/to/6F_adjusted.1Mb_bins.binCov_matrix.bed.gz
```

The output from this ploidy estimation script is a large number of files and plots, all of which summarize and/or visualize per-sample ploidy per chromosome.  

#### Optional: Calcuating per-sample and binwise coverage distribution properties  
Median values can be evaluated with ```medianCoverage.R```. This tool operates in two modes (1) medians per bin across all samples, and (2) medians per sample across all bins. For each, median coverages are returned when considering all data and also when excluding bins (or samples) with zero coverage. See [documentation of medianCoverage.R](https://github.com/RCollins13/WGD#mediancoverager)  for details of this process.  

--- 

### CNV visualization and annotation  
#### Companion visualization tool: CNView  
[CNView](https://github.com/RCollins13/CNView) is a companion tool for WGD that can visualize, score, and annotate CNVs directly from ___raw___ ```makeMatrix.sh``` output. More details on CNView can be found [here](http://biorxiv.org/content/early/2016/04/20/049536). If you use CNView, please cite [Collins et al., 2016](http://biorxiv.org/content/early/2016/04/20/049536).

---  
## Script Documentation (incomplete; work in progress)  
---  
### binCov.py
Iterates through a single chromosome of a bam file and calculates either nucleotide or physical coverage in regularly segmented bins.
```
usage: binCov.py [-h] [-S] [-C] [-I INDEX_PATH] [-z] [-n NORM_OUT]
                 [-b BINSIZE] [-m {nucleotide,physical}] [-x BLACKLIST]
                 [-v OVERLAP] [--oldBT]
                 bam chr cov_out

Calculates non-duplicate primary-aligned binned coverage of a chromosome from
an input BAM file

positional arguments:
  bam                   Input bam (or sam/cram, but requires appropriate flag)
  chr                   Contig to evaluate
  cov_out               Output bed file of raw coverage

optional arguments:
  -h, --help            show this help message and exit
  -S, --SAM             Input file is in sam format
  -C, --CRAM            Input file is in cram format
  -I INDEX_PATH, --index_path INDEX_PATH
                        Bam/cram index file
  -z, --gzip            Gzip output files bed files
  -n NORM_OUT, --norm_out NORM_OUT
                        Output normalized coverage
  -b BINSIZE, --binsize BINSIZE
                        Bin size, in bp (default: 1000)
  -m {nucleotide,physical}, --mode {nucleotide,physical}
                        Evaluate nucleotide or physical coverage (default:
                        nucleotide)
  -x BLACKLIST, --blacklist BLACKLIST
                        BED file of regions to ignore
  -v OVERLAP, --overlap OVERLAP
                        Maximum tolerated blacklist overlap before excluding
                        bin
  --oldBT               Flag to indicate if you are using a BEDTools version
                        pre-2.24.0
```  
**Usage Notes:**  
- Input bam file must be coordinate-sorted and indexed.  
- Only non-duplicate primary-aligned reads or proper pairs are considered for 'nucleotide' and 'physical' mode, respectively.  
- Normalized coverage is raw coverage per bin divided by median of all non-zero, non-blacklisted bins on the same contig.  
- Bins will be ignored automatically if they share at least ```-v``` percent overlap by size with blacklisted regions (```-x``` or ```--blacklist```).  
- Currently uses ```bedtools coverage``` syntax assuming ```bedtools``` version 2.24.0 or later (i.e. ```-a``` is bins and ```-b``` is reads; this was reversed starting in ```bedtools v2.24.0```). Specify option ```--oldBT``` to revert to bedtools syntax pre-2.24.0.  

---  

### WG_binCov.sh
Wrapper for serialized execution of binCov.py across multiple chromosomes for an individual sample  
```
usage: WG_binCov.sh [-h] [-b BINSIZE] [-m MODE] 
                    [-L CONTIGS] [-x BLACKLIST] [-v OVERLAP] 
                    BAM ID OUTDIR

Wrapper for serialized execution of binCov.py across multiple chromosomes

Positional arguments:
  BAM     Input bam
  ID      Sample ID
  OUTDIR  Output directory

Optional arguments:
  -h  HELP         Show this help message and exit
  -b  BINSIZE      Bin size in bp (default: 1000)
  -m  MODE         Evaluate physical or nucleotide coverage (default: nucleotide)
  -L  CONTIGS      List of contigs to evaluate (default: all contigs in bam header)
  -x  BLACKLIST    BED file of regions to ignore
  -v  OVERLAP      Maximum tolerated blacklist overlap before excluding bin
```  
**Usage Notes:**  
- Contents of arguments are not checked; it's up to the user to ensure they're feeding appropriate files and values.  

---  

### makeMatrix.sh
Helper tool to automate creation of sorted coverage matrices from binCov.py output bed files. Wraps ```bedtools unionbedg```. Note that, in most cases, multiple chromosome outputs from binCov.py from the same sample should be concatenated before being passed to this tool.  
```
usage: makeMatrix.sh [-h] [-z] [-o OUTFILE] SAMPLES

Helper tool to automate creation of sorted coverage matrices from
binCov.py output bed files

Positional arguments:
  SAMPLES     List of samples and coverage files (tab-delimmed)

Optional arguments:
  -h  HELP      Show this help message and exit
  -z  GZIP      Gzip output file
  -o  OUTFILE   Output file (default: stdout)
```
**Usage Notes:**  
- Does not sanity check to ensure that sample IDs are identical between matrices, or that samples were included in the same order when passed to makeMatrix.sh.

---  

### medianCoverage.R
Helper tool to calculate median coverages per-bin across all samples or per-sample across all bins. Medians are reported both with and without considering zero-coverage bins or samples. Uses a sorted coverage matrix as can be generated by makeMatrix.sh (see above).
```
usage: medianCoverage.R [-h] [-b/--binwise] covMatrix.bed OUTFILE

Helper tool to calcualte √arious medians from binCov matrices

Positional arguments:
  covMatrix.bed   path to binCov coverage matrix
  OUTFILE         full path to desired output file

Optional arguments:
  -h  HELP      Show this help message and exit
  -b/--binwise  Compute medians of all samples per bin [default: median of all bins per sample]
```  
**Usage Notes:**  
- NOTE: Automatically downsamples to 1M bins (in per-sample mode) or 500 samples (in per-bin mode) to increase computational efficiency. Will add an option at a later date to disable downsampling.
- DEV NOTE: will be extended to calcuate other measurements, such as standard deviation, mean, median absolute deviation, quartiles, and min/max

---  

### compressCov.sh
Helper tool to compress bincov output (individual sample files or mutli-sample matrices) into larger bin sizes. By default, median values per new bin are reported.
```
usage: compressCov.sh [-h] [-z] [-n] [-o OUTFILE] INPUT RATIO

Helper tool to automate compression of raw binCov.py output bed files or
bed-style coverage matrices into larger bin sizes

Positional arguments:
  INPUT     path to binCov.py bed file or bed-stype matrix
  RATIO     compression ratio

Optional arguments:
  -h  HELP        Show this help message and exit
  -z  GZIP        Gzip output file
  -s  SUM         Report sum (default: report median)
  -o  OUTFILE     Output file (default: stdout)
```  
**Usage Notes:**  
- Compression ratio (RATIO) must be a positive integer. New bins will be automatically instantiated with a size equal to RATIO times the current binsize. New bins that do not overlap any previous bins (e.g. bins were excluded due to blacklisting (-x) during binCov.py) will not be reported.
- Input file does not have to be split by chromosome.
