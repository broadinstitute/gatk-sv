---
title: Quotas
description: Guidance on quota limits and increases in GCP
sidebar_position: 5
slug: quotas
---

When running GATK-SV in Google Cloud, users may need to be aware of quotas. 
Quotas restrict resource usage in each Google Cloud project. Reaching quota 
limits can slow down workflow execution or result in workflow failures. 
Generally, default quotas are appropriate for processing up to a few thousand 
samples with GATK-SV, but quota increases may be required for larger cohorts. 

### Monitoring quotas and requesting increases

For information on how to monitor quota usage and request quota increases, refer 
to [this article](https://support.terra.bio/hc/en-us/articles/360029071251-How-to-troubleshoot-and-fix-stalled-workflows).

### Relevant quotas for GATK-SV

The following quotas may need to be increased in order to efficiently process
large cohorts with GATK-SV:

- Regional managed instance groups
- In-use regional external IPv4 addresses
- CPUs
- Persistent Disk Standard (GB)
- VM Instances
- Preemptible CPUs
- Queries per minute per region
- Operation read requests per minute per region
- Google Egress Bandwidth per second per region
- Persistent Disk SSD (GB)


#### Regional managed instance groups

The largest bottleneck initially is likely to be regional managed instance groups (RMIGs).
The default quota for RMIGs in `us-central1` is 1,250. Cromwell operates one VM per RMIG, 
so this is equivalent to maximum of 1,250 concurrently-running VMs.

#### IP addresses

Large quota increases for external IP addresses (in-use regional external IPv4 addresses) 
are difficult to obtain. Therefore, GATK-SV workflows requiring the highest degree of 
parallelization (such as `GatherSampleEvidence`, when applied to a large cohort) may 
use internal IP addresses instead. A task will use an internal IP address if the task's 
`runtime` block in the WDL contains `noAddress: true`. 

Tasks using internal IP addresses do not exhaust the external IP address quota. However,
they are limited by the size of the VPC network. By default, this is 4,096 IPs. 
This limit is not listed on the quotas page in the cloud console; it is on the VPC networks 
page, but Terra users do not have access to this page. Terra support can increase the size
of the VPC network upon request. When the limit is reached for internal IP addresses, 
error messages will mention `IP space of <project's subnetwork> is exhausted. Insufficient free IP addresses`.
