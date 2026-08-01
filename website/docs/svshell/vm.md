---
title: VM Setup
description: Setting up the VM for SV-Shell
sidebar_position: 1
slug: vm-setup
---

import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';

SV-Shell is a platform-agnostic implementation of the single-sample pipeline: 
it runs end-to-end in a single Docker image, 
so it works on any platform that can run Docker images, 
such as Azure, AWS, GCP, and on-premises infrastructure.



## Running SV-Shell

The following steps walk through provisioning a VM, 
installing dependencies, 
and running the single-sample pipeline end-to-end inside the SV-Shell Docker image.


1. Create a Linux VM.
The VM should have at least 30 GB RAM and 50GB disk space on the OS disk and ~500GB mounted disk space.
We recommend `Standard L2aos v4` type on Azure, and `c4-standard-4-lssd` on GCP.

2. Connect to the VM

3. Install the needed packages on the VM.

    <Tabs
      groupId="cloud"
      defaultValue="gcp"
      values={[
        { label: 'GCP', value: 'gcp' },
        { label: 'Azure', value: 'azure' },
      ]}>
      <TabItem value="gcp">
        ```bash
        sudo apt-get update && sudo apt-get install -y tmux
        ```
      </TabItem>
      <TabItem value="azure">
        ```bash
        sudo apt-get update && sudo apt-get install -y tmux && \
        wget https://aka.ms/downloadazcopy-v10-linux && \
        tar -xvf downloadazcopy-v10-linux && \
        sudo cp ./azcopy_linux_amd64_*/azcopy /usr/bin/
        ```
      </TabItem>
    </Tabs>

4. Configure mounted disks
We assume the disks are mounted to the following path:

    ```shell
    # /mnt/disks/gatk-sv/sv-shell/
    ```

    <Tabs
      groupId="cloud"
      defaultValue="gcp"
      values={[
        { label: 'GCP', value: 'gcp' },
        { label: 'Azure', value: 'azure' },
      ]}>
      <TabItem value="gcp">
        You may use the following docs on mounting and partitioning disks or merging multiple SSDs into a single logical partition:
        https://docs.cloud.google.com/compute/docs/disks/add-local-ssd#formatmultiple
      </TabItem>
      <TabItem value="azure">
        You may use the script we have developed for this purpose:
        ```shell
        # list nvme:
        lsblk -o NAME,SIZE,TYPE,FSTYPE,MOUNTPOINT,MODEL
    
        # then use this script to format and mount them
        wget https://gist.githubusercontent.com/VJalili/9f3bbbb9f34a41ae6a2099d25f75008f/raw/98724fddc62aedefc77f58183fdd80bbf11745b0/format-single-nvme.sh .
    
        bash format-single-nvme.sh <nvme ID>
        sudo chmod a+w /mnt/disks/gatk-sv/sv-shell
        ```
    
        Or manually configure depending on the VM family type.
        For L-series VMs that need RAID, you may take the following steps.
    
        ```shell
        # first use this command to get the list of nvmes:
        lsblk -o NAME,SIZE,TYPE,FSTYPE,MOUNTPOINT,MODEL
    
        # next partition them:
        sudo apt-get update && sudo apt-get install -y mdadm && \
        yes | sudo mdadm --create /dev/md0 --level=0 --name=nvme_raid --raid-devices=3 /dev/nvme1n1 /dev/nvme2n1 /dev/nvme3n1 && \
        sudo mkfs.ext4 -F /dev/md0 && \
        sudo mkdir -p /mnt/disks/gatk-sv/sv-shell/ && \
        sudo mount /dev/md0 /mnt/disks/gatk-sv/sv-shell/ && \
        sudo chown -R $USER:$USER /mnt/disks/gatk-sv/sv-shell/
        ```
      </TabItem>
    </Tabs>

    Ensure the disk is mounted/configured correctly.

    ```shell
    df -h
    
    # expected output
    Filesystem       Size  Used Avail Use% Mounted on
    /dev/root         24G  2.6G   21G  11% /
    ...
    /dev/md0         737G  2.1M  700G   1% /mnt/disks/gatk-sv/sv-shell
    ```

5. Localize files

    ```shell
    mkdir -p /mnt/disks/gatk-sv/sv-shell/inputs && cd /mnt/disks/gatk-sv/sv-shell/inputs
    ```

    <Tabs
      groupId="cloud"
      defaultValue="gcp"
      values={[
        { label: 'GCP', value: 'gcp' },
        { label: 'Azure', value: 'azure' },
      ]}>
      <TabItem value="gcp">
        ```shell
        # Remove tracker files. This ensures that if you re-run the command after interrupting
        # a previous run, gcloud starts from a fresh tracker; otherwise you may see errors like:
        # ⠹ERROR: Expecting value: line 1 column 1 (char 0) 2.8MiB/s
        rm -rf ~/.config/gcloud/surface_data/storage/tracker_files/*
    
        time gcloud storage cp [the bucket containing all the data] .
        ```
      </TabItem>
      <TabItem value="azure">
        ```shell
        azcopy copy '[the blob containing data with the SAS token appended]' '/mnt/disks/gatk-sv/sv-shell/inputs/' --recursive
        ```
      </TabItem>
    </Tabs>

6. Install Docker on the VM:

    ```shell
    sudo apt update && \
    sudo apt install -y apt-transport-https ca-certificates curl software-properties-common && \
    curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo gpg --dearmor -o /usr/share/keyrings/docker-archive-keyring.gpg && \
    echo \
        "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/docker-archive-keyring.gpg] https://download.docker.com/linux/ubuntu \
        $(lsb_release -cs) stable" | sudo tee /etc/apt/sources.list.d/docker.list > /dev/null && \
    sudo apt update && \
    sudo apt install -y docker-ce docker-ce-cli containerd.io && \
    sudo usermod -aG docker $USER && newgrp docker
    ```

7. Run the Docker image:

    ```shell
    docker run --rm -t -d --entrypoint /bin/bash --platform linux/amd64 \
      -v "/mnt/disks/gatk-sv/sv-shell/inputs/:/inputs/" \
      -v "/mnt/disks/gatk-sv/sv-shell/wd:/wd/" \
      us.gcr.io/broad-dsde-methods/vjalili/sv-shell:52aeeb52
    ```

8. Start and enter a tmux session

    ```shell
    tmux
    ```

9. Get the ID of the running Docker container using `docker ps`

10. Exec into the container:

    ```shell
    docker exec -it {CONTAINER_ID} /bin/bash
    ```

11. Config the single-sample pipeline.

    ```shell
    cd /opt/sv_shell
    export SV_SHELL_BASE_DIR="/wd"
    export TMPDIR="/wd/tmp"
    mkdir -p /wd/tmp
    ```

12. To collect stats:

    ```shell
    export STATS_DIR="${SV_SHELL_BASE_DIR}/stats_$(date +'%Y%m%d_%H%M%S')"
    mkdir -p "${STATS_DIR}"
    /usr/lib/sysstat/sadc 1 "${STATS_DIR}/" &
    SADC_PID=$!
    ```

13. Start the single-sample pipeline

    ```shell
    { time bash single_sample_pipeline.sh sample_inputs/single_sample_pipeline.json; } \
      > >(tee "${SV_SHELL_BASE_DIR}/stdout.log") \
      2> >(tee "${SV_SHELL_BASE_DIR}/stderr.log" >&2)
    ```

14. If you were running stats, you can extract them to a TSV file as follows:

    ```shell
    kill ${SADC_PID}
    sleep 2
    ls -tr "${STATS_DIR}"/sa?? | xargs -I {} sar -A -f {} > "${STATS_DIR}/stats.txt"
    ```
