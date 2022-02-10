This repository includes scripts to run minGQ on gnomad samples.

The original minGQ script has now been splitted into three parts:

*Module07MinGQPart1CollectData.wdl*: 

	- the first script to run for minGQ; 
	- collects information from cleanvcf to build minGQ model; 
	- this script can parallelize across chromosomes; 
	- json files are kept under **example_jsons/** with **minGQ.gnomad_v3.Step1.** in file name.

*Module07MinGQPart2TrainModel.wdl*: 

	- the second script to run for minGQ; 
	- integrates data across the genome to build global minGQ model;  
	- outputs across all chromosomes from Part1 are integrated to build the minGQ model in this script;
	- json files are kept under **example_jsons/** with **minGQ.gnomad_v3.Step2.** in file name.

*Module07MinGQPart3ApplyModel.wdl*: 

	- last script to run for minGQ; 
	- applies the trained model from Part2 to filter SVs; 
	- this wdl can parallelize across chromosomes;

