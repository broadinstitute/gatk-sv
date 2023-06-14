---
title: Manual Deployment
description: Build and Publish Images
sidebar_position: 3
---

If you are contributing to the GATK-SV codebase, specifically focusing on 
enhancing tools, configuring dependencies in Dockerfiles, or modifying GATK-SV scripts 
within the Docker images, it is important to build and test the Docker images locally. 
This ensures that the images are successfully built and function as intended. 
Additionally, if you wish to host the images in your own container registry, 
you will need to follow these steps. 
To simplify the build process, we have developed a Python script 
that automates the image building, and publishing to your container registry. 
This section provides a detailed guide on using this script.

