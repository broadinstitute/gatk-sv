#!/usr/bin/env python

import argparse
import sys
import os
import os.path
import time
import tempfile
import shutil
import pprint
# noinspection PyPackageRequirements
from termcolor import colored
import subprocess


this_script_path = os.path.dirname(os.path.abspath(__file__))


class ProjectBuilder:
    """
    class to track dependencies, control build and push of entire job
    """
    github_org = 'broadinstitute'
    github_repo = 'gatk-sv'
    # mapping from target to its dependencies
    #   each dependency is either None, or a mapping from each dependency name to the docker ARG it is passed via
    #   currently each image has zero or one dependencies, but multiple dependencies are allowed
    dependencies = {
        'delly': None, 'manta': None, 'melt': None, 'wham': None, 'sv-base-mini': None,
        'samtools-cloud': {'sv-base-mini': "MINIBASE_IMAGE"}, 'sv-base': {'sv-base-mini': "MINIBASE_IMAGE"},
        'cnmops': {'sv-base': "SVBASE_IMAGE"},
        'sv-pipeline-base': {'sv-base': "SVBASE_IMAGE"},
        'sv-pipeline': {'sv-pipeline-base': "SV_PIPELINE_BASE_IMAGE"},
        'sv-pipeline-children-r': {'sv-pipeline-base': "SV_PIPELINE_BASE_IMAGE"},
        'sv-pipeline-hail': {'sv-pipeline': "SV_PIPELINE_IMAGE"},
        'sv-pipeline-updates': {'sv-pipeline': "SV_PIPELINE_IMAGE"},
        'sv-pipeline-rdtest': {'sv-pipeline-children-r': "SV_PIPELINE_BASE_R_IMAGE"},
        'sv-pipeline-qc': {'sv-pipeline-children-r': "SV_PIPELINE_BASE_R_IMAGE"}
    }
    temporary_images = frozenset({})  # currently no temporary images are needed, but they can be used
    non_public_images = frozenset({'melt'})
    images_built_by_all = frozenset(dependencies.keys()).difference({"melt"})
    accepted_target_values = frozenset(dependencies.keys()).union({"all"})

    def __init__(self, project_arguments, launch_script_path=this_script_path):
        # Todo: we can also check auth to push to GCR
        if project_arguments.gcr_project is not None:
            os.system("docker login")

        self.project_arguments = project_arguments
        self.launch_script_path = launch_script_path
        self.working_dir = None
        self.build_priority = {}

    def get_build_priority(self, target_name):
        if target_name not in self.build_priority:
            build_deps = ProjectBuilder.dependencies[target_name]
            self.build_priority[target_name] = 0 if build_deps is None \
                else 1 + max((self.get_build_priority(build_dep) for build_dep in build_deps.keys()), default=0)
        return self.build_priority[target_name]

    def get_ordered_build_chain_list(self):
        targets_to_traverse = set(self.project_arguments.targets)
        targets_to_build = set()
        while targets_to_traverse:
            target = targets_to_traverse.pop()
            if self.project_arguments.no_force_rebuild and \
                    ImageBuilder.image_is_built(f"{target}:{self.project_arguments.image_tag}"):
                continue  # this target is already built, so it and its dependencies and be ignored
            targets_to_build.add(target)
            dependencies = ProjectBuilder.dependencies[target]
            if dependencies is not None:
                new_targets_to_traverse = set(dependencies.keys()).difference(targets_to_build)
                if self.project_arguments.skip_base_image_build:  # only build temporary dockers
                    new_targets_to_traverse.intersection_update(ProjectBuilder.temporary_images)
                targets_to_traverse.update(new_targets_to_traverse)

        build_chain = sorted(targets_to_build, key=self.get_build_priority)
        return tuple(build_chain)  # immutable, at least attempt

    def go_to_workdir(self):
        tmp_dir_path = None
        # local mode
        if self.project_arguments.staging_dir is None:
            # working_dir is grandparent directory of launch script path
            self.working_dir = os.path.dirname(os.path.dirname(self.launch_script_path))
        else:
            # if staging is required, mkdir, cd, and pull
            tmp_dir_path = tempfile.mkdtemp(prefix=self.project_arguments.staging_dir).rstrip('/') + '/'
            connect_mode = "git@github.com:" if self.project_arguments.use_ssh else "https://github.com"
            clone_target = connect_mode + "/" + ProjectBuilder.github_org + "/" + ProjectBuilder.github_repo + ".git"
            if os.system(f"git clone {clone_target} {tmp_dir_path}") != 0:
                raise RuntimeError(f"Failed to clone {clone_target}.")
            self.working_dir = tmp_dir_path
        os.chdir(self.working_dir)

        # checkout desired hash or tag, if building remotely
        if self.project_arguments.remote_git_tag is not None:
            git_checkout_cmd = "git checkout tags/" + self.project_arguments.remote_git_tag
            if os.system(git_checkout_cmd) != 0:
                raise ValueError(f"The provided git tag [{self.project_arguments.remote_git_tag}] does not exist")
        elif self.project_arguments.remote_git_hash is not None:
            git_checkout_cmd = "git checkout      " + \
                self.project_arguments.remote_git_hash
            if os.system(git_checkout_cmd) != 0:
                raise ValueError(f"The provided git hash [{self.project_arguments.remote_git_hash}] does not exist")

        print("Working directory: " + os.getcwd())
        return tmp_dir_path

    @property
    def remote_docker_repos(self):
        return () if self.project_arguments.gcr_project is None \
            else (f"us.gcr.io/{self.project_arguments.gcr_project}",)

    @property
    def built_images(self):
        # note: need to separate out format part of string because docker format conflicts with f-string
        return get_command_output(
            f"docker images *:{self.project_arguments.image_tag} " + "--format '{{.Repository}}:{{.Tag}}'"
        ).split()

    def build(self):
        if "all" in self.project_arguments.targets:
            self.project_arguments.targets = ProjectBuilder.images_built_by_all
        expanded_build_targets = self.get_ordered_build_chain_list()
        print("Building the following targets in order:")
        print(expanded_build_targets)
        print(colored('#################################################', 'magenta'))

        for target_name in expanded_build_targets:
            a = colored("Building image ", "grey")
            b = colored(target_name + ":" + self.project_arguments.image_tag, "yellow", attrs=['bold'])
            c = colored(" ...", "grey")
            print(a, b, c)

            dependencies = ProjectBuilder.dependencies[target_name]
            build_time_args = {} if dependencies is None else {
                arg: f"{image_name}:{self.project_arguments.image_tag}" for image_name, arg in dependencies.items()
            }

            ImageBuilder(target_name, self).build(build_time_args, self.project_arguments.no_force_rebuild)
            print(colored('#################################################', 'magenta'))

    def push(self):
        for image in self.built_images:
            ImageBuilder(image, self).push()

    def build_and_push(self):
        print("Project args:")
        pprint.pprint(vars(self.project_arguments))
        print("")
        possible_tmp_dir_path = self.go_to_workdir()

        try:
            # start docker daemon, if one hasn't been started yet
            os.system("open --background -a Docker && while ! docker system info > /dev/null 2>&1; do sleep 1; done")
            self.build()
            self.push()
            print(colored('BUILD PROCESS SUCCESS!', 'green'))
        finally:
            if not self.project_arguments.skip_cleanup:
                self.cleanup(possible_tmp_dir_path)

    def cleanup(self, tmp_dir_path):
        # clean dangling images (i.e. those "<none>" images), stopped containers, etc
        os.system("docker system prune -f")

        # clean intermediate image that are not to be pushed
        for local_image in self.built_images:
            if ImageBuilder(local_image, self).name in ProjectBuilder.temporary_images:
                os.system(f"docker rmi --force {local_image}")

        # "rm -rf" staging dir, if was specified
        if (self.project_arguments.staging_dir is not None) and (tmp_dir_path is not None):
            os.chdir(self.launch_script_path)  # first cd back to launch_script_path
            shutil.rmtree(tmp_dir_path)


class ImageBuilder:  # class for building and pushing a single image
    def __init__(self, name, project_builder):
        if ':' in name:
            self.name, self.tag = name.split(':', 1)
        else:
            self.name = name
            self.tag = project_builder.project_arguments.image_tag
        self.project_builder = project_builder

    @property
    def working_dir(self):
        return self.project_builder.working_dir

    @property
    def local_image(self):
        return f"{self.name}:{self.tag}"

    @property
    def remote_docker_repos(self):
        return self.project_builder.remote_docker_repos

    @staticmethod
    def image_is_built(local_image):
        images = {
            image for image in get_command_output(
                'docker images --format "{{.Repository}}:{{.Tag}}"'
            ).split('\n')
        }
        return local_image in images

    def build(self, built_time_args_dict, no_force_rebuild=False):
        if no_force_rebuild and ImageBuilder.image_is_built(self.local_image):
            print(f"skipping build of {self.local_image} because it is already built")
            return
        # standard build command
        docker_build_command = "docker build --progress plain \\\n    "
        docker_build_command += "-f " + f"{self.working_dir}/dockerfiles/{self.name}/Dockerfile" + " \\\n    "
        docker_build_command += "--tag " + self.local_image + " \\\n    "
        # parse extra args list
        for key, value in built_time_args_dict.items():
            docker_build_command += "--build-arg " + key + "=" + value + " \\\n    "

        will_push = 0 != len(self.remote_docker_repos) and any(
            e is not None for e in self.remote_docker_repos)
        docker_build_command += "--squash . " if will_push else ". "

        # build and time it
        print(docker_build_command)
        start_time = time.time()
        if os.system(docker_build_command) != 0:
            raise RuntimeError(f"Failed to build image {self.local_image}")
        elapsed_time = time.time() - start_time
        elapsed_min, elapsed_sec = divmod(elapsed_time, 60)
        print("Time spent on docker build:")
        print(f"{elapsed_min} minutes, {elapsed_sec} seconds")

    def push(self):
        if self.name in ProjectBuilder.temporary_images:
            return  # don't push temporary images
        self._push(self.tag)
        if self.project_builder.project_arguments.update_latest:
            self._push("latest")

    def _push(self, remote_tag):
        for rep in self.remote_docker_repos:
            # do not push images with very restrictive licenses
            if (self.name in ProjectBuilder.non_public_images) and (not rep.startswith('us.gcr.io')):
                print(colored(f"Refusing to push non-public image {self.name} to {rep}", "red"))
                continue

            remote_image = f"{rep}/{self.name}:{remote_tag}"
            docker_tag_command = f"docker tag {self.local_image} {remote_image}"
            docker_push_command = f"docker push {remote_image}"
            print(docker_tag_command)
            if os.system(docker_tag_command) != 0:
                raise RuntimeError(f"Failed to tag image ({remote_image}) for pushing to remote")
            print(docker_push_command)
            if os.system(docker_push_command) != 0:
                raise RuntimeError(f"Failed to push image {remote_image}")


def get_command_output(command):
    """
    Execute shell command. Raise exception if unsuccessful, otherwise return string with output
    """
    sub_p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    with sub_p.stdout as pipeIn, sub_p.stderr as pipeErr:
        results = pipeIn.read().decode('utf-8')
        err = pipeErr.read().decode('utf-8')
    if err:
        raise RuntimeError('Error executing %s:\n%s' % (command, err[:-1]))
    return results


def __parse_arguments(args_list):
    parser = argparse.ArgumentParser(
        description='=' * 50 + "\nBuilding docker images for GATK-SV pipeline v1.\n" + '=' * 50,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    # required arguments
    required_args_group = parser.add_argument_group('Required', 'required arguments')
    required_args_group.add_argument(
        '--targets', nargs='+', type=str, required=True,
        help='the sub project docker(s) you want to build (note "all" does not include melt)'
    )
    required_args_group.add_argument('--image-tag', type=str, required=True,
                                     help='tag to be applied to all images being built')
    # to build from local or remote git tag/hash values
    git_args_group = parser.add_argument_group('Mutex args', 'remote git tag/hash values (mutually exclusive)')
    git_mutex_args_group = git_args_group.add_mutually_exclusive_group()

    git_mutex_args_group.add_argument(
        '--remote-git-tag', type=str,
        help='release tag on Github; this indicates pulling from Github to a staging dir'
    )
    git_mutex_args_group.add_argument(
        '--remote-git-hash', type=str,
        help='a hash value on Github; this indicates pulling from Github to a staging dir'
    )
    # to build from remote/Github (staging required)
    remote_git_args_group = parser.add_argument_group('Remote git',
                                                      'args involved when building from remote git tags/hashes')

    remote_git_args_group.add_argument(
        '--staging-dir', type=str, help='a temporary staging directory to store builds; required only when pulling '
                                        'from Github, ignored otherwise'
    )
    remote_git_args_group.add_argument('--use-ssh', action='store_true', help='use SSH to pull from github')
    # flag to turn on push to Dockerhub and/or GCR
    docker_remote_args_group = parser.add_argument_group(
        'Docker push', 'controlling behavior related pushing dockers to remote repos'
    )
    docker_remote_args_group.add_argument('--gcr-project', type=str,
                                          help='GCR billing project to push the images to. If not given, the '
                                               'built image(s) will not be push to GCR.')
    docker_remote_args_group.add_argument('--update-latest', action='store_true',
                                          help='also update \"latest\" tag in remote docker repo(s)')
    # flag to turn off git protection (default mode is refusing to build when there are untracked files and/or
    # uncommitted changes)
    parser.add_argument(
        '--disable-git-protect', action='store_true',
        help='disable git check/protect when building from local files (will use uncommited changes to build)'
    )
    parser.add_argument('--skip-base-image-build', action='store_true',
                        help='skip rebuild of the target\'s base image(s). Assumes that the base image(s) already '
                             'exist with same tag.')
    parser.add_argument('--no-force-rebuild', action='store_true',
                        help='Do not rebuild docker images if the exact image and tag already exist.')
    parser.add_argument('--skip-cleanup', action='store_true',
                        help='skip cleanup after successful and unsuccessful build attempts.')

    # parse and consistency check
    if len(args_list) <= 1:  # no arguments, print help and exit with success
        parser.parse_args(['-h'])
        sys.exit(0)
    parsed_args = parser.parse_args(args_list[1:])

    # if passed targets are in the accepted values
    for tar in parsed_args.targets:
        if tar not in ProjectBuilder.accepted_target_values:
            raise ValueError("\"" + tar + "\" not in allowed target values")

    if "all" in parsed_args.targets:
        if 1 != len(parsed_args.targets):
            raise ValueError("when \"all\" is provided, no other target values allowed")

    # if "use_ssh" flag is turned on, remote git tag/hash should be provided
    if parsed_args.use_ssh is True:
        if (parsed_args.remote_git_tag is None) and (parsed_args.remote_git_hash is None):
            raise ValueError("\"use_ssh\" is specified but remote git tag/hash is not")

    # if remote git tag/hash and/or is specified, staging dir should be specified
    if (parsed_args.remote_git_tag is not None) or (parsed_args.remote_git_hash is not None):
        if parsed_args.staging_dir is None:
            raise ValueError("remote git tag/hash is specified but staging_dir is not")

    # if requesting to update "latest" tag in remote docker repo(s), remote git release tag must be specified
    if parsed_args.update_latest is True:
        if parsed_args.remote_git_tag is None:
            raise ValueError("publishing \"latest\" docker images requires a remote Github release tag")

    # if there are uncommitted changes when building from local files, raise exception
    if parsed_args.staging_dir is None and not parsed_args.disable_git_protect:
        s = os.popen("git status -s | wc -l | tr -d ' ' | tr -d '\n'").read()
        ret = int(s)
        if 0 != ret:
            raise ValueError(
                "Current directory has uncommitted changes or untracked files. Cautiously refusing to proceed."
            )
    return parsed_args


def main(arguments):
    project_args = __parse_arguments(arguments)
    ProjectBuilder(project_args).build_and_push()


if __name__ == "__main__":
    main(sys.argv)
