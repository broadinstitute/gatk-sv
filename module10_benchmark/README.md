This repository includes scripts that annotate cleanvcf from gnomadV3 with the features that are necessary for boost filtering. Function of each script is explained as follow:

**SplitPerSampleGTGQ.wdl** - script designed to split GT and GQ information of each sample from cleanvcf. XXX samples are recommended to be processed together. An example of the input json file can be found at: **

**AnnotateILPBFeaturesBed.wdl** -  script designed to annotate HGSV samples that have PacBio seqeunces, output will annotate each sample with basic features, including: GT and GQ splitted from vcf, overlap by raw algorithms (manta, wham, melt, depth), overlap with matched PacBio and Bionano SVs, support from PacBio sequences as evaluated by VaPoR. An example of the input json file can be found at: *./example_jsons/AnnotateILPBFeaturesBed.PB_samples.cleanvcf.json*

**AnnotateILFeaturesBed.wdl** - script designed to annotate all gnoamd samples with basic features required by boost model, including GT and GQ splitted from vcf and overlap by raw algorithms (manta, wham, melt, depth). A maximum of XXX samples can be processed as a single batch. An example of the input json file can be found at: **


