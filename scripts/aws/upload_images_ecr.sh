#!/bin/bash
set -e -x
# This script will download the images required for the pipeline from GCP/Dockerhub and will upload them to ECR.
# This is needed to avoid timeouts/limits on the GCP/docker side.
# This also helps fetch the images quickly on AWS while running the jobs.
# Make sure docker services are running on the machine from where this is being run.

function usage() {
  printf "Usage: \n \
    %s -r <AWS_REGION> -e <ECR_REPO_NAME> \n \
    <AWS_REGION> \t the AWS region where the images need to be uploaded. Ex: us-east-1, us-east-2 etc. \n \
    <ECR_REPO_NAME> \t The ECR repository name where all the images will be uploaded." "$1"
}

if [[ "$#" == 0 ]]; then
  usage "$0"; exit 0;
fi

#################################################
# Parsing arguments
#################################################
while getopts "r:e:" option; do
  case "$option" in
    r) AWS_REGION="$OPTARG" ;;
    e) ECR_REPO_NAME="$OPTARG" ;;
    *) usage "$0" && exit 1 ;;
  esac
done

if [ -z "$AWS_REGION" ] ; then
    echo "Checking for AWS Region"
    usage "$0"
    exit 1
fi

if [ -z "$ECR_REPO_NAME" ] ; then
    echo "Checking for ECR Repo Name"
    usage "$0"
    exit 1
fi



# This function will check the exit code and return status.
function check_status_of_command(){
    if [ $1 -ne 0 ]
    then
        echo "The process has failed at $2. Please check and resubmit."
		exit 11;
    else
        echo "Completed : $2 ."
    fi
}

# Download the json file locally
# This can be replaced by 
# cp ../../inputs/values/dockers.json docker.json
wget https://github.com/broadinstitute/gatk-sv/blob/master/inputs/values/dockers.json\?raw\=true -O docker.json
check_status_of_command $? "Downloading https://github.com/broadinstitute/gatk-sv/blob/master/input_values/dockers.json Locally"


#Get Account ID and login to ECR
aws configure set region ${AWS_REGION}
export AWS_ACCOUNT_ID=$(aws sts get-caller-identity --output text --query 'Account')
aws ecr get-login-password --region ${AWS_REGION} | docker login --username AWS --password-stdin ${AWS_ACCOUNT_ID}.dkr.ecr.${AWS_REGION}.amazonaws.com
check_status_of_command $? "ECR Login for Account : ${AWS_ACCOUNT_ID}"

# Create ECR repo if not exists
aws ecr create-repository --repository-name ${ECR_REPO_NAME} || true

# Pull Image from GCR/Docker Hub and Push Image to ECR
pull_and_push_image()
{
    docker pull $1
    check_status_of_command $? "Docker Pull : ${1}"
    docker image tag $1 ${AWS_ACCOUNT_ID}.dkr.ecr.${AWS_REGION}.amazonaws.com/${ECR_REPO_NAME}:${2}
    check_status_of_command $? "Docker tagging : ${1}"
    docker push ${AWS_ACCOUNT_ID}.dkr.ecr.${AWS_REGION}.amazonaws.com/${ECR_REPO_NAME}:${2}
    check_status_of_command $? "Docker Push : ${AWS_ACCOUNT_ID}.dkr.ecr.${AWS_REGION}.amazonaws.com/${ECR_REPO_NAME}:${2}"
}

while read line
do
    # Ignoring google cloud sdk image and lnux dockerhub as it will be downloaded on runtime.
    if [[ "$line" == *"_docker"* ]] && [[ "$line" != *"cloud_sdk_docker"* ]] && [[ "$line" != *"linux_docker"* ]] && [[ "$line" != *"melt_docker"* ]]; then
        echo "Running for : $line"
        ecr_image_tag=`echo $line | cut -d '"' -f2`
        gcr_image=`echo $line | cut -d '"' -f4`
        echo "The GCR Image : $gcr_image"
        echo "ECR Image Tag : $ecr_image_tag"
        pull_and_push_image $gcr_image $ecr_image_tag
    fi
done < docker.json

echo "The process is complete."