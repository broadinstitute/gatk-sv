# Update the below variables
S3_OR_FSX=$1
BROAD_REF_PATH="${S3_OR_FSX}/reference/broad-ref"
CRAM_BAM_FILE_PATH="${S3_OR_FSX}/bams"
HAPLOTYPE_GVCF_PATH="${S3_OR_FSX}/reference/gvcf"
GATK_SV_RESOURCES_PATH="${S3_OR_FSX}/reference/gatk-sv-resources"
BATCH_DEF_PATH="${S3_OR_FSX}/reference/batch_sv.test_large.qc_definitions.tsv"
AWS_ACCOUNT_ID=$(aws sts get-caller-identity --output text --query 'Account')
AWS_REGION=`aws configure get region`
ECR_REPO_NAME="sv-pipeline"

# Install docker
sudo yum install -y jq
sudo amazon-linux-extras install -y docker
sudo usermod -a -G docker ec2-user
sudo service docker start
sudo chkconfig docker on
sudo curl -L https://github.com/docker/compose/releases/download/1.22.0/docker-compose-$(uname -s)-$(uname -m) -o /usr/local/bin/docker-compose
sudo curl -L https://github.com/docker/compose/releases/latest/download/docker-compose-$(uname -s)-$(uname -m) -o /usr/local/bin/docker-compose
sudo chmod +x /usr/local/bin/docker-compose
docker-compose version
sudo chmod 777 /var/run/docker.sock

# Download the gatk-sv github repo
# This will be hard-coded to a particular release/tag if Broad is unable to maintain it.
# git clone https://github.com/spatel-gfb/gatk-sv.git

# Create the required code and reference files.
cd gatk-sv
mkdir gatk_run 
python3 -m pip install jinja2
BASE_DIR=$(pwd)
GATK_SV_ROOT=$(pwd)
CLOUD_ENV="aws.gatk_sv"
echo '{ "google_project_id": "broad-dsde-methods", "terra_billing_project_id": "broad-dsde-methods" }' > inputs/values/${CLOUD_ENV}.json
python3 scripts/inputs/build_inputs.py ${BASE_DIR}/inputs/values ${BASE_DIR}/inputs/templates/test ${BASE_DIR}/inputs/build/ref_panel_1kg/test  -a '{ "test_batch" : "ref_panel_1kg", "cloud_env" : "'$CLOUD_ENV'" }'
cp inputs/build/ref_panel_1kg/test/GATKSVPipelineBatch/GATKSVPipelineBatch.json
cp inputs/build/ref_panel_1kg/test/GATKSVPipelineBatch/GATKSVPipelineBatch.json gatk_run/

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

# Update and Copy the aws config file.
cd ../../scripts/aws/
# Update the aws_GATKSVPipelineBatch.ref_panel_1kg.json json with the correct values as per variables defined
array=( BROAD_REF_PATH CRAM_BAM_FILE_PATH HAPLOTYPE_GVCF_PATH GATK_SV_RESOURCES_PATH BATCH_DEF_PATH AWS_ACCOUNT_ID AWS_REGION ECR_REPO_NAME )
for i_var in "${array[@]}"
do
	i_val=${!i_var}
    # Below might need -i '' if running on mac.
	sed -i "s/${i_var}/${i_val}/g" aws_GATKSVPipelineBatch.json
done
cp opts.json aws_GATKSVPipelineBatch.json ../../gatk_run/

# Upload the images to ECR
sh upload_images_ecr.sh -r ${AWS_REGION} -e ${ECR_REPO_NAME}


echo "IMPORTANT : Kindly compare the below :
    - BROAD : gatk-sv/gatk_run/GATKSVPipelineBatch.json
    - AWS : gatk-sv/gatk_run/aws_GATKSVPipelineBatch.json

And update the missing params in AWS json from BROAD json with correct AWS paths/account id/aws region and at the same location as of BROADs.
"