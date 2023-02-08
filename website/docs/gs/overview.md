---
title: Overview
description: Setup and options for running GATK-SV.
sidebar_position: 0
---

# Overview

GATK-SV is a highly-scalable cloud-native pipeline for structural variation discovery 
on Illumina short-read whole-genome sequencing (WGS) data.
The pipeline genotypes structural variations using Docker-based tools, modularized in 
different components, and orchestrated using Cromwell.

The pipeline runs in two modes: [Cohort](/docs/execmodes/cohort) and [single-sample](/docs/execmodes/singlesample). 
