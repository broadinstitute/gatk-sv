Thank you for your interest in contributing to GATK-SV!

When it comes to open-source, every contribution makes the 
software better for everyone and is appreciated by the 
community. To express our gratitude for your contribution, 
we would like to provide you with easy-to-follow steps to get started.

# Code of Conduct

In the interest of fostering an open and welcoming environment, we as 
contributors and maintainers pledge to make participation in our 
project and our community a harassment-free experience for everyone, 
regardless of age, body size, disability, ethnicity, sex characteristics, 
gender identity and expression, level of experience, education, 
socio-economic status, nationality, personal appearance, race, religion, 
or sexual identity and orientation. 

GATK-SV and everyone participating in it is governed by the 
[Broad Institute Code of Conduct](https://github.com/broadinstitute/.github/blob/develop/CODE_OF_CONDUCT.md).
By participating, you are expected to uphold this code.


# What Should I Know Before I Get Started?

GATK-SV is a _cloud-native pipeline_ written in [Workflow Description Language](https://openwdl.org/) (WDL) 
that orchestrates joint genotyping structural variations using open-source tools and custom scripts. 
The tools and scripts are distributed in Docker containers, and the pipeline is executed using a 
[Cromwell](https://cromwell.readthedocs.io/en/stable/) 
server running on commercial cloud platforms. You may interface with a Cromwell server using 
[Cromshell](https://github.com/broadinstitute/cromshell) or [Terra platform](https://terra.bio). 
We currently publish Docker images on Google Container Registry (GCR) and Azure Container Registry (work-in-progress).

To study your data using the GATK-SV pipeline, please 
refer to the [documentation](README.md) for details.

# I have a question! 

Please refer to the [documentation](README.md). 
If you still have questions, please start a new issue under the 
[Issues](https://github.com/broadinstitute/gatk-sv/issues) tab. 


# How Can I Contribute?

## Reporting Bugs

We appreciate your bug reports; bugs are tracked as 
[GitHub issues](https://github.com/broadinstitute/gatk-sv/issues). 
Before submitting a bug report, please search the issues tab, 
as you may find a similar issue reported. Please submit a 
bug report if you do not find a similar issue.

### How do I submit a bug report?

The following guidelines help submit bug reports that are 
easier to understand and reproduce by the maintainers and 
the community and help find related reports.

- Use a clear and descriptive title for the issue; 
- Provide a minimally reproducible example; 
- Provide specific details on how the observed and expected behaviors differ.


## Your First Code Contribution

Thank you for considering code contributions to GATK-SV. 
All code contributions to the GATK-SV repository are made through 
[pull requests (PR)](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-pull-requests), 
and you may take the following steps to prepare and submit your PRs.

### Fork GATK-SV (one time)
A [fork](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/about-forks) 
is your copy of a GitHub repository that you can make your changes 
without affecting the source.

You may create a fork of [broadinstitute/gatk-sv](https://github.com/broadinstitute/gatk-sv) 
following [these instructions](https://docs.github.com/en/get-started/quickstart/fork-a-repo#forking-a-repository).

### Clone GATK-SV fork (one time)

Once you have forked [broadinstitute/gatk-sv](https://github.com/broadinstitute/gatk-sv), 
you need to download your fork to your computer so you can implement any changes. 
This process is commonly known as `clone`ing a GitHub repository. You may follow 
[these instructions](https://docs.github.com/en/get-started/quickstart/fork-a-repo#cloning-your-forked-repository) 
on cloning your fork of GATK-SV, or in a nutshell, you may run the following commands.

* Clone your fork
    ```shell
    git clone https://github.com/YOUR_GITHUB_USERNAME/gatk-sv .
    cd gatk-sv
    ```

* Add a reference to `broadinstitute/gatk-sv` so you can keep your fork in sync in the future.

    ```shell
    git remote add upstream https://github.com/broadinstitute/gatk-sv
    ```

### Create a feature branch (once per feature)

In a loose sense, a [_branch_](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-branches) 
is an isolated copy of the code in your forked repository, where the changes 
you implement in one branch does not affect the code on the other branches. 
It is a common practice to create one branch per each feature you would like 
to implement, and do not implement features on the `main` branch of the repository. 
It is mainly because you can use the `main` branch to synchronize your fork with 
the upstream (i.e., `broadinstitute/gatk-sv`) and get a fresh copy of the up-to-date codebase to start your feature development. 

You may create a branch using [GitHub](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-and-deleting-branches-within-your-repository) 
or [git](https://git-scm.com/book/en/v2/Git-Branching-Basic-Branching-and-Merging). 
In a nutshell, you may create a branch as the following. 

* Ensure you are on the `main` branch.  
    ```shell
    git checkout main
    ```

* Create a branch.
    ```shell
    git checkout -b YOUR-BRANCH-NAME
    ```

* Synchronize your branch with `broadinstitute/gatk-sv`. 
    ```shell
    git fetch upstream
    git merge upstream/main
    ```

* Push your branch to GitHub.
    ```shell
    git push –set-upstream origin YOUR-BRANCH-NAME
    ```

### Implement a new feature or fix a bug

You may open the code base in any text editor or integrated development environment (IDE). 
We recommend using [PyCharm](https://www.jetbrains.com/pycharm/) or 
[Visual Studio](https://visualstudio.microsoft.com/vs/community/), 
which are both freely available. Please follow the following guidelines 
when making changes to the codebase.

- **Coding Style**: GATK-SV is composed of code written in various programming languages, 
including WDL, Python, Bash, and R. Please follow the styling guidelines of each programming 
language; e.g., [PEP-8](https://peps.python.org/pep-0008/) for Python.

- **Documentation**: Please add concise and descriptive documentation to your changes.



### Commit changes
_Committing_ changes means tracking them via _git_. You may submit a change in 
one commit or split it into multiple related commits. Though there is no hard 
rule on splitting your commits, there are 
[best practices](https://github.blog/2022-06-30-write-better-commits-build-better-projects/) 
that we highly recommend. In GATK-SV, we 
[`squash and merge`](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/incorporating-changes-from-a-pull-request/about-pull-request-merges#squash-and-merge-your-commits) 
pull requests, which squashes all the commits into a single commit and 
merge on the `main` branch.

You may commit changes to git using its command line interface (CLI) from a terminal 
or use the graphical user interface (GUI) of the IDE you are using. 
You may use the following command if you are using a terminal 
(details [here](https://www.atlassian.com/git/tutorials/saving-changes/git-commit), 
or follow the guidelines on making commits in 
[PyCharm](https://www.jetbrains.com/help/pycharm/commit-and-push-changes.html#commit) 
and [Visual Studio](https://learn.microsoft.com/en-us/visualstudio/version-control/git-make-commit?view=vs-2022)).

```shell
git commit -m "commit message"
```

### Make a pull request

The commits you created are tracked on the feature branch of your fork. 
To bring these changes to GATK-SV, you need to create a pull request. 

- You may follow these guidelines on 
[creating a pull request](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request?tool=webui&platform=windows).

- Choose a concise and descriptive title for your pull request 
  that briefly summarizes your proposed improvements. 
  Avoid unnecessary details such as your initials or branch names in the title.

- Provide concise and descriptive information on the expected behavior.

- Thoroughly test your changes and provide minimally reproducible steps to test your code.

After your pull request is reviewed and approved, your commits can be squash-committed 
to the `main` branch of GATK-SV.


# License
By contributing your code to the GATK-SV GitHub repository, you agree 
to license your contribution under the [BSD 3-Clause license](LICENSE.TXT).
