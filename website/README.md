# Website

The GATK-SV [technical documentation website](https://broadinstitute.github.io/gatk-sv/)
is built using [Docusaurus](https://docusaurus.io/). The website is built
from the Markdown (`.md`) files available in the `website/docs/` directory.
The website is built into the [GATK-SV CI/CD](https://github.com/broadinstitute/gatk-sv/blob/main/.github/workflows/docs.yml),
hence, any changes to the `website/` is automatically tested, built, and deployed.


## Installation
You may build and run it for local development using the following commands.

### Prerequisite

- Install [Node.js](https://nodejs.org/en/download/);
- Install [Git LFS](https://docs.github.com/en/repositories/working-with-files/managing-large-files/installing-git-large-file-storage)
since GATK-SV tracks images using `git lfs`. After you `git clone` the GATK-SV repository, run the 
following command to fetch the `LFS`-tracked files. 
    ```
    $ git lfs install
    $ git lfs fetch --all
    $ git lfs pull
    ```
- Mac M1 users: You may need to [set IPv6 configuration](https://support.apple.com/guide/mac-help/use-ipv6-on-mac-mchlp2499/mac)
to `Link-Local Only` for `npm` to successfully install packages.  

### Install 

```
$ cd website/
$ npm install
```

This command installs necessary packages and creates `package-lock.json`
for the reproducibility of deployment.

### Build

```
$ npm run build
```

This command generates static content into the `build` directory and can be served using any static contents hosting service.

### Start

```
$ npm run start
```

This command starts a local development server and opens up a browser window. Most changes are reflected live without having to restart the server.
