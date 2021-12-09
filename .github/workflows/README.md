# CI/CD

The continuous integration (CI) and continuous delivery (CD) pipeline of 
GATK-SV is developed on [GitHub Actions](https://docs.github.com/en/actions/learn-github-actions/understanding-github-actions).
The CI/CD pipeline is defined via multiple _workflows_ where each is 
a `.yml` file stored under the `.github/workflows` directory. The workflows
are triggered automatically when a pull request (PR) is issued or merged.
The workflows automate testing, building, and deploying the pipeline, 
and they currently cover the following areas. 

- Lint Python scripts (`pytest.yaml`): 
using [flake8](https://pypi.org/project/flake8/) asserts if the Python scripts
follow the PEP-8 style guides;
- Test, build, and publish docker images using
[`build_docker.py`](/scripts/docker/build_docker.py)
(`sv_pipeline_docker.yml`).

## Test, Build, and Publish Docker images

The GATK-SV Docker images are built and published using
the [`build_docker.py`](/scripts/docker/build_docker.py), which
is documented at [this README](/scripts/docker/README.md) and can be
executed locally.

The [`Docker Images workflow`](sv_pipeline_docker.yml) (DIW) automates the
test, build, and publication of GATK-SV Docker images using `build_docker.py`,
such that, the images are built when a PR is issued against the repository,
and published to Google Cloud Platform (GCP) Container Registry 
(GCR) (at `us.gcr.io`) when the PR is merged.

### Workflow Layout

The DIW consists of three [_jobs_](https://docs.github.com/en/actions/learn-github-actions/workflow-syntax-for-github-actions#jobs):
1. `Determine Build Args`.
This job determines the arguments to be used by the `build_docker.py` script,
specifically:
   - Given the size and the number of GATK-SV Docker images, DIW builds and 
     publishes only the Docker images affected by the changes introduced in 
     a PR. Accordingly, first the files changed between the `HEAD` and 
     the `BASE` commits of the PR'd branch are determined using `git diff`
     (for details, please refer to the in-line comments for 
     the step `Determine Commit SHAs` in [`DIW`](sv_pipeline_docker.yml)), 
     and then the affected images are determined. These images are used 
     as the values of `--targets` argument of the `build_docker.py` script.
   - A step to compose a tag for the Docker images in the `DATE-HEAD_SHA_8`
     template, where `DATE` is `YYYYMMDD` extracted from the time stamp 
     of the last commit on the PR'd branch, and `HEAD_SHA_8` is the first 
     eight characters of its commit SHA. For instance `20211201-86fe06fd`.


2. `Test Images Build`. This job is triggered when a commit
  is pushed to the PR'd branch; it builds docker images determined by
  the `Determine Build Args` job. This job fails if the building of 
  the Docker images was unsuccessful. The Docker images built by this job
  will not be published to GCR and are discarded as the job succeeds.


3. `Publish`. This job is triggered when a PR is merged to a commit 
  is pushed to the `master` branch. Similar to the `Test Images Build` job,
  this job builds Docker images and fails if the build process was 
  unsuccessful. However, in addition, this job pushes the built images
  to GCR. To authorize access to GCR, this job assumes a GCP service 
  account (SA) with read and write access to the GCR registry. The secrets 
  related to the SA are defined as 
  [encrypted environment secrets](docs.github.com/en/actions/security-guides/encrypted-secrets).


### Setup the `Deploy` Environment
_This section describes configuring the `Deploy` environment to be used
by the `Publish` job and is intended for the edification of repository admins.

An SA is used to authorize DIW to access GCR. (A future extension may 
adopt an [OpenID Connect [OIDC]-based](https://docs.github.com/en/actions/deployment/security-hardening-your-deployments/about-security-hardening-with-openid-connect)
authentication and authorization). In order to assume the SA, the `Publish`
job needs the SA secrets (e.g., `private key` and `client email`) and 
`project name`. This information is defined in a [GitHub environment](https://docs.github.com/en/actions/deployment/targeting-different-environments/using-environments-for-deployment),
and is exposed to the `Publish` job as encrypted environment secrets.
The encrypted secrets are decrypted in the environment context and are 
not exposed to the user code (if the norms of best practices are followed).
GitHub's environment secrets are a subset of repository-wide secrets, 
which is a subset of organization-level secrets. We encrypt SA credentials 
as GitHub's environment secrets as they allow pausing the execution of any 
action that accesses the environment until it is approved by assigned 
individuals.

In order to set up the `Deploy` environment, you may take the following steps:

1. [Create an SA on GCP IAM](https://cloud.google.com/iam/docs/creating-managing-service-accounts#creating).
   For simplicity, you may assign the service account the `Editor` role.
   However, in order to follow the principles of the least privilege, 
   you may assign the `Storage Object Admin`, `Storage Legacy Bucket Writer`,
   and `Storage Object Viewer` as the minimum required permissions
   ([ref](https://cloud.google.com/container-registry/docs/access-control)).


3. Get the service account's keys by going to the 
   [`Service Accounts page`](https://console.cloud.google.com/iam-admin/serviceaccounts)
   and selecting the above-created service account and going to the `KEYS` tab.
   Then click on the `ADD KEY` button, and choose `Create new key`. In the 
   pop-up window, select `JSON` type and click on the `CREATE` button. It
   will download a JSON file containing the secrets required to assume the 
   service account.


4. Base64 encode the service account's secrets in the JSON format as the 
   following.

   ```shell
   openssl base64 -in service-account.json -out service-account.txt
   ```

5. Create an environment following [these steps](https://docs.github.com/en/actions/deployment/targeting-different-environments/using-environments-for-deployment#creating-an-environment)
   *and name the environment `Deploy`*. 


6. Create the following two encrypted secrets in the `Deploy` environment 
   using [these steps](https://docs.github.com/en/actions/security-guides/encrypted-secrets#creating-encrypted-secrets-for-an-environment):
   - `name`: `GCP_PROJECT_ID`; `value`: the ID of the GCP project 
     under which you will use the GCR registry.
   - `name`: `GCP_GCR_SA_KEY`; `value`: the above-created base64 encoding 
     of the SA's secrets. After you set this encrypted secret, we 
     recommend that you delete both the `.json` and `.txt` files 
     containing SA's secrets.


7. [Optional] Under the `Environment protection rules` on the `Deploy` 
   environment's configuration page, you may check the `Required reviewers`
   checkbox and assign maintainers who can approve the execution of the 
   instances of the jobs that require access to the `Deploy` environment. 


### Review pending deployments

Once the `Deploy` environment is set up, and the `Required reviewers`
option under the section `Environment protection rules` is checked, 
with every push to the `master` branch (e.g., merging a PR), the
DIW execution will pause at the `Publish` job with the following 
message: 

> Waiting for review: Deploy needs approval to start deploying changes.

If enabled, any `Required reviewers` will see the following
additional link that they can click to approve or reject running the 
job.

> Review pending deployments
