---
title: Runtime attributes
description: Runtime attributes
sidebar_position: 6
---

GATK-SV is implemented as a set of WDLs designed to run on Google Cloud Platform. Computations are broken up into 
a set of tasks that are carried out in a particular order on cloud virtual machines (VMs), each of which 
possesses a limited set of resources. These include the following primary components:

- CPU cores
- Random access memory (RAM)
- Disk storage

GATK-SV attempts to request the appropriate amount of resources for each task and in some cases employs mathematical
models to predict optimal requirements that minimize cost. However, if the actual computations performed require more
than the requested resources on the VM, the task may fail.

:::info
Most tasks in GATK-SV are tuned for cohorts with ~150 samples. Most tasks should scale automatically, but running larger 
cohorts may lead to slowdowns or errors due to insufficient resource allocation.
:::

In addition to VM resources, Terra and Cromwell both support automatic retries of tasks in case of ephemeral errors, 
as well as requesting [preemptible VMs](https://cloud.google.com/compute/docs/instances/preemptible) at a significant discount.

### Setting runtime attributes

GATK-SV exposes optional parameters for manually tuning VM resource allocation, automatic retries, and the use of preemptible 
instances. These parameters are all a custom struct called `RuntimeAttr`, which is defined in 
[Structs.wdl](https://github.com/broadinstitute/gatk-sv/blob/main/wdl/Structs.wdl) as:

```
struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}
```

Users encountering errors due to insufficient VM resources or wishing to adjust automatic and preemptible retries may 
modify the corresponding parameters. Users should inspect the WDL source to determine whether a given `RuntimeAttr` 
parameter corresponds to a specific task.

### Further reading

For more information on identifying errors and setting runtime parameters, please refer to the 
[Troubleshooting FAQ](/docs/troubleshooting/faq).
