## Genotype filtering optimization  

_Note: this directory contains instructions to run the revised version of the workflow previously known as minGQ, now known as Module07FilterGTs_ 

The original minGQ workflow has now been divided into three sequential steps:

1. **../wdl/Module07FilterGTsPart1CollectData.wdl**:  
```
- Collects data for all trios present in an input VCF to train a genotype filtering model in Step 2 
- This step can be parallelized across chromosomes
- Input json files are kept under example_jsons/  with "minGQ.gnomad_v3.Step1." in file name [note: these reflect XZ's original implementation]
```

2. **../wdl/Module07FilterGTsPart2TrainModel.wdl**: 
```
- Integrates trio data across all chromosomes output by Step 1 to build global genotype filtering model  
- json files are kept under example_jsons/ with "minGQ.gnomad_v3.Step2." in file name [note: these reflect XZ's original implementation]
```

3. **../wdl/Module07FilterGTsPart3ApplyModel.wdl**: 
```
- Filters genotypes using the model trained in Step 2
- This step can be parallelized across chromosomes
- Recommended to apply directly to the sharded, PCR-split VCFs output by step 1 (for reasonable runtimes)
- json files are kept under example_jsons/ with "minGQ.gnomad_v3.Step3." in file name [note: these reflect XZ's original implementation]
```

---  

Developers note: this codebase is a work-in-progress and there may be some inconsistencies between XZ's original implementation and RLC's revised implementation  
