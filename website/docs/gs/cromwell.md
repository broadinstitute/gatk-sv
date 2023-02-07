---
title: Cromwell
description: Run GATK-SV on a Cromwell server.
sidebar_position: 6
---

[Cromwell](https://github.com/broadinstitute/cromwell) is a workflow management system
that takes a workflow (e.g., a workflow written in [Workflow Description Language (WDL)](https://openwdl.org)), 
its dependencies and input data, and runs it on a given platform 
(e.g., 
[Azure](https://cromwell.readthedocs.io/en/stable/backends/Azure/), 
[GCP](https://cromwell.readthedocs.io/en/stable/backends/Google/), or 
[AWS](https://cromwell.readthedocs.io/en/stable/backends/AWSBatch/)). 
In order to run a workflow on Cromwell, you need a running instance of 
Cromwell that is available in two forms: [Server and stand-alone mode](https://cromwell.readthedocs.io/en/stable/Modes/).

In general, you may use a managed Cromwell server maintained by your 
institute or host a self-managed server, or run Cromwell as a standalone Java application.
The former is ideal for large scale execution in a managed environment, 
and the latter is useful for small scale and isolated WDL development.

:::info
Due to its dependency on cloud-hosted resources and large-scale execution needs,
we currently do not support running the entire GATK-SV pipeline using 
Cromwell as a [stand-alone Java application](https://cromwell.readthedocs.io/en/stable/Modes/#run) 
:::

# Cromwell Server

There are two option to communicate with a running Cromwell server: 
[REST API](https://cromwell.readthedocs.io/en/stable/tutorials/ServerMode/), and
[Cromshell](https://github.com/broadinstitute/cromshell) which is a command line tool
to interface with a Cromwell server. We recommend using Cromshell due to its simplicity 
of use. This documentation is explained using Cromshell, but the same steps can be 
taken using the REST API.

## Setup Cromwell

1. Enable the compute engine for your project by navigating to the
  `“Compute Engine”` tab of the cloud console (making sure your new project is selected). 
  It will take a few minutes to load. To check that Compute Engine has been configured, 
  navigate to `“IAM & Admin” > “Service accounts”`. You should see a line for  
  `<project-number>-compute@developer.gserviceaccount.com`.

2. Add the Cromwell service account for the relevant Cromwell server to your project 
  as an Editor via the IAM & Admin page of the cloud console. 
  The cromwell service account needs project-level access to the project you are 
  running compute on in order to launch VMs through the compute engine. 
  In this case bucket access is not sufficient.

3. Set up a cromwell execution bucket. Options:

   - Create a new bucket in your new project and set a time-to-live. 
     Then make sure to specify the bucket path in your cromwell config 
     JSON with the “jes_gcs_root” option as it is not the default bucket for the cromwell server.

   - If you want to use an existing bucket in a separate project 
     (e.g., the default cromwell execution bucket), grant read/write permissions for 
     the bucket to the compute engine service account for your new project 
     (`<project-number>-compute@developer.gserviceaccount.com`). 
     The role Storage Admin should suffice.


4. Grant appropriate permissions to the Compute Engine default service account 
  (for your new project) **and** Cromwell service account. For buckets that they 
  may write to, use the role Storage Object Admin (read + write access); 
  for buckets for which only read-access is necessary, use the role Storage Object Viewer. 
  If granting permissions on a project level, Editor is an appropriate role for 
  service accounts. Both service accounts will need appropriate permissions 
  for the following buckets/objects:

   - Cromwell execution bucket;

   - All input data (buckets or objects). If input data is in Terra but 
     you want to use cromshell to process it (rather than Terra), 
     see this section on how to add service accounts to Terra

   - Final output bucket (if using final_workflow_outputs_dir workflow option);

   - Monitoring script (if using monitoring_script);

   - Docker images in Google Container Repository (most GATK-SV 
     images are public, in which case no action is needed, but Melt is not public).

5. To allow the use of the PAPIv2 genomics API by cromwell, 
  go to the Genomics API page (even though it says “deprecated,” 
  this is the right page), select your new project in the dropdown 
  on the top menu bar, and click “Enable API.” There is no need to 
  add any credentials for this API to use it with the default compute engine service account.

6. In order to bill compute engine costs to the right project, 
  you must specify the “google_project” option in your cromwell config JSON. 
  Specify the project ID (all lowercase) as opposed to the project name.
  Format: `“google_project”: “project-id”`.

7. To allow the creation of compute engine VMs, create a network, 
  subnet, and firewall rule limiting access to Broad devices 
  (or devices on the Broad VPN) by following these instructions 
  but with the following caveats:

   - The network name broad-allow is used as an example throughout 
     the instructions. When creating the network, subnetwork, 
     and firewall rule, specify your own network name rather 
     than broad-allow. For instance, project-name-net.

   - Skip the subnet;

   - When creating the firewall rule: `--allow=tcp`, skipp `--target-tags` argument.

   - Close Cloud Shell after you are done by clicking the X on the right side

   - Another resource on this topic: https://dsp-security.broadinstitute.org/cloud-security/google-cloud-platform/securing-the-network 


8. Go to https://console.cloud.google.com/iam-admin/labels and select 
  your project from the dropdown menu, then add two labels with keys 
  “wdl-network-name” and “wdl-subnetwork-name” and the value is the 
  network name (the one you just created!) for both.

9. When submitting your workflow, add "backend": "PAPIv2-vpc" to your cromwell workflow options config JSON.


## Setup Cromshell

1. Install Cromshell:

   ```shell
   $ brew install broadinstitute/dsp/cromshell 
   ```
  
2. Configure Cromshell to communicate with a Cromwell server, and make sure 
  it can communicate with the server. For instance, if the server is behind a 
  firewall (i.e., accessible to only the members of your institute), you may 
  need to be connected to a VPN.

  The first time you use cromshell, it will ask you to set your default cromwell server URL. 
  If you want to use a different server going forward, or temporarily, update 
  the cromwell server you’re using by either:

   - Running the command `export CROMWELL_URL=https://<server.org>` 
     (replace with relevant server URL), to change the cromwell server 
     for the rest of your bash session.

   - Adding the above command to your .bash_profile or .bashrc file 
     (and running source ~/.bashrc) to set it as the default server 
     for every new bash session.

   - To set the cromwell server for just one command at a time, 
     add `CROMWELL_URL=https://<server.org>` (or the desired URL) 
     to the beginning of your cromshell command, like the below example:

     ```shell
     $ CROMWELL_URL=https://<server.org> cromshell -t60 metadata [workflow-id] > meta.json
     ```



3. Create a Cromwell config JSON file. An example basic config file is 
  shown below (more options here and here):

   ```json
   {
     "workflow_failure_mode": "ContinueWhilePossible",
     "monitoring_script": "gs://gatk-sv-resources/resources/cromwell_monitoring_script2.sh",
     "read_from_cache" : true,
     "write_to_cache" : true
   }
   ```
   
   - The monitoring script has to be in Google Cloud when running a 
     remote Cromwell workflow. If you don’t have access to the copy of the 
     script in the example above, you can upload your own to another Google 
     Cloud bucket. The script is here in the gatk-sv git repo.
   
   - Other options that are particularly useful include:
      - `“jes_gcs_root”: “gs://bucket-name”`: specify the location of the cromwell execution bucket;
      - `“final_workflow_outputs_dir”: “gs://bucket-name”`: specify a bucket to which cromwell 
        should copy the top-level workflow outputs (important if the execution bucket has a 
        TTL but you want to retain the outputs);
      - `"google_project": "talkowski-sv-gnomad"`: specify the google billing project. 
        Important if you want to bill a project (for compute costs) other than the default 
        one for the cromwell server. Some finagling of permissions may be required.

4. Make a directory (outside the git repo, or just don’t push it to the remote repository) 
  for testing and cd into that directory.

5. Copy in your WDL of interest, all the dependent WDLs (or just all the WDLs, because 
  that’s easier), and the input JSON. Let’s say you want to run the test_small test set 
  through ModuleX.wdl, which is in the /wdl/ directory of the git repo - then your 
  commands might look something like this: 

   ```shell
   $ cp /path/to/repo/gatk-sv/test/moduleX/ModuleX.test_small.json .
   $ cp -r /path/to/repo/gatk-sv/wdl/ .
   ```

6. Create a zip file of all the dependent WDLs.

   ```shell
   $ zip -j dependencies.zip /path/to/repo/gatk-sv/wdl/*.wdl
   ```
   
7. Submit the job using Cromshell, specifying the workflow you want to run, its inputs, 
  the Cromshell configuration, and the dependencies.


   ```shell
   $ cromshell submit ModuleX.wdl ModuleX.test_small.json /path/to/cromwell_config.json dependencies.zip
   ```
   
   - When you submit the job, Cromshell will run the womtool validation script on it first 
     (if womtool is in your path - see docs), to validate the primary WDL and input JSON. 
     If your submission fails that validation, it will print an error message, 
     and you can fix the errors and try again.

   - If the job is submitted successfully, it will display a turtle and a job ID 
     (a long hash of numbers and lowercase letters). 

8. Some useful commands for monitoring your job include (see docs for more).

   ```shell
   # See the status
   $ cromshell status <job ID>
   
   # View a timing diagram in browser to track workflow progress, You can use workflow sub-IDs for this
   $ cromshell timing <job ID>
   
   # Get metadata (after job completes)
   $ cromshell -t 60 metadata <job ID>
   
   # Get a stripped-down version 
   $ cromshell slim-metadata <job ID>
   ```