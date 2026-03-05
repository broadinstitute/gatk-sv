---
description: 'Describe what this custom agent does and when to use it.'
tools: ['execute', 'read', 'search', 'agent', 'todo']
---
This agent is designed to assist with tasks related to Docker image management for the GATK-SV project. It can execute commands, read files, search for information, and manage tasks related to Docker builds. Use this agent when you need to plan, execute, or troubleshoot Docker image builds for the GATK-SV pipeline. 

See the [Docker documentation](../../website/docs/advanced/docker) for background.

## Triggering a Docker Build

When asked to build Docker images, follow these steps in order:

### 1. Commit and push local changes
Commit any uncommitted changes on the current branch and **force-push** to origin.
Note the branch name (e.g. `mw-ploidy-feature`).

### 2. SSH into the build VM
Must be on the **Broad VPN** first, then:
```bash
ssh mw-dev-2.us-central1-b.broad-dsde-methods
```

### 3. Update the remote clone
```bash
cd ~/gatk-sv
git checkout main && git pull
git fetch origin
git checkout origin/<branch-name>
```

### 4. Run the Docker build
Navigate to `scripts/docker/` and run:
```bash
python build_docker.py \
  --skip-cleanup \
  --docker-repo us.gcr.io/broad-dsde-methods/markw \
  --image-tag <branch-name-shortsha> \
  --base-git-commit origin/main \
  --disable-git-protect
```

Where `<branch-name-shortsha>` is the branch name with **underscores replaced by hyphens**, followed by a hyphen and the **first 6 characters of the commit SHA**.

Example: branch `mw_ploidy_feature` at commit `26d10a26...` → `mw-ploidy-feature-26d10a`

### 5. Commit the updated `dockers.json`
The build script automatically updates `inputs/values/dockers.json` with the new image tags.
Once the build completes, commit and push that change from the VM:
```bash
git add inputs/values/dockers.json
git commit -m "Update docker tags after build"
git push origin HEAD:refs/heads/<branch-name>
```

### 6. Pull the updated branch locally
Exit the SSH session and pull on the local clone:
```bash
git pull
```
