# Update the below variables
BROAD_REF_PATH="S3 or FSx mount path where broad-ref reference files exist"
CRAM_BAM_FILE_PATH="S3 or FSx mount path where Sample BAM or CRAM (wth index) files exist"
HAPLOTYPE_GVCF_PATH="S3 or FSx mount path where Haplotype caller gvcf files exist"
GATK_SV_RESOURCES_PATH="S3 or FSx mount path where gatk-sv-resources files exist"
BATCH_DEF_PATH="S3 or FSx mount path where batch_sv.test_large.qc_definitions.tsv or equivalent file resides"
AWS_ACCOUNT_ID="AWS Account Id"
AWS_REGION="AWS Region in which operating"
ECR_REPO_NAME="ECR Repo Name needed ex: 'sv-pipeline'"


# Download the gatk-sv github repo
# This will be hard-coded to a particular release/tag if Broad is unable to maintain it.
git clone https://github.com/broadinstitute/gatk-sv.git

# Create the required code and reference files.
cd gatk-sv
mkdir gatk_run
cd gatk_run
cp -r ../wdl .
cd wdl; 
# MAKE MELT FLAGS FALSE. Search for "use_melt" and change flag to false. 
# `find . -type f | xargs grep "Boolean use_melt"` can be used to search files.
sed -i "s/Boolean use_melt = true/Boolean use_melt = false/g" GATKSVPipelineBatch.wdl
sed -i "s/Boolean use_melt = true/Boolean use_melt = false/g" GATKSVPipelineSingleSample.wdl
# Update TinyResolve CPU/MEM for AWS in order to run it on newer instance as we have seen Container Pull issues while it runs with other jobs.
# Increasing the CPU/MEM to 16 will ensure a new Batch EC2 is spun up coz rest other parallel jobs are using below 8 cpus.
sed -i "s/cpu_cores: 1/cpu_cores: 16/g;s/mem_gb: 3.75/mem_gb: 16/g" TinyResolve.wdl
zip dep.zip *.wdl
cd ../../scripts/inputs/ 
python3 -m pip install jinja2   # Needed in AWS as EC2 doesnt have it installed.
python3 build_inputs.py ../../input_values ../../input_templates/GATKSVPipelineBatch.ref_panel_1kg.json.tmpl . -a '{ "ref_panel" : "ref_panel_1kg_v2", "test_batch" : "test_batch_large"}'
cp GATKSVPipelineBatch.ref_panel_1kg.json ../../gatk_run/

# Update and Copy the aws config file.
cd ../scripts/aws/
# Update the aws_GATKSVPipelineBatch.ref_panel_1kg.json json with the correct values as per variables defined
array=( BROAD_REF_PATH CRAM_BAM_FILE_PATH HAPLOTYPE_GVCF_PATH GATK_SV_RESOURCES_PATH BATCH_DEF_PATH AWS_ACCOUNT_ID AWS_REGION ECR_REPO_NAME )
for i_var in "${array[@]}"
do
	i_val=${!i_var}
    # Below might need -i '' if running on mac.
	sed -i "s/${i_var}/${i_val}/g" aws_GATKSVPipelineBatch.ref_panel_1kg.json
done
cp opts.json aws_GATKSVPipelineBatch.ref_panel_1kg.json ../../gatk_run

# Upload the images to ECR
sh upload_images_ecr.sh -r ${AWS_REGION} -e ${ECR_REPO_NAME}


echo "IMPORTANT : Kindly compare the below :
    - BROAD : gatk-sv/gatk_run/GATKSVPipelineBatch.ref_panel_1kg.json 
    - AWS : gatk-sv/gatk_run/aws_156_samples_GATKSVPipelineBatch.ref_panel_1kg.json 

And update the missing params in AWS json from BROAD json with correct AWS paths/account id/aws region and at the same location as of BROADs.
"
