#!/usr/bin/env python

import argparse
import sys
import os
import fnmatch
import json
import time
import tempfile
import shutil
import pprint
# noinspection PyPackageRequirements
from termcolor import colored
import subprocess
from typing import Union, Iterable, Mapping, Optional, List, Tuple, Set, Dict
from types import MappingProxyType

program_operation_description = """
The high-level view of how build_docker.py works:
    1) Primary target selection. This is the images we specifically
       want to re-build. Determined either:
        a) manually (specified via --targets), or
        b) automatically (via changes detected from git, specified
           via --base-git-commit and --current-git-commit)
    2) Determine which images must be built. Done by finding images
       that depend (directly or indirectly on targets).
       If --no-force-rebuild is set and a target or dependent image
       is already built, exclude it from being rebuilt, and stop the
       chain of dependencies (presumably this is most useful as part
       of a debugging process where several Dockerfiles are being
       tweaked).
    3) Load current default docker images from dockers json specified
       by --input-json
    4) For each image that must be built:
        a) Actually build the image. If it depends on other images,
           supply the dependent image tag as specified by the input
           json.
        b) If there is a remote image repo specified, push the newly
           built image to that repo and update the dockers json with
           the remote image.
        c) If there is no remote image repo, update the dockers json
           with the local image.
    5) If a remote image repo is specified and there are pre-existing
       local images in the dockers json that were not built this run,
       push them now and update the dockers json with their remote
       name.
    6) Write out the updated dockers json to the file specified by
       --output-json

When you run the program normally, it prints statements describing its
actions at each step, including printing out changes to dockers.json.
When you run it with --dry-run it does all the same things, except it:
    a) never builds an image
    b) never pushes an image
    c) doesn't actually alter the output dockers json
Instead it prints the changes that it would make if you ran without
--dry-run.
"""


class Paths:
    this_script_folder = os.path.dirname(os.path.abspath(__file__))
    gatk_sv_path = os.path.dirname(os.path.dirname(this_script_folder))
    dockers_json_path = os.path.join(gatk_sv_path, "inputs", "values", "dockers.json")
    dev_null = "/dev/null"


class ImageDependencies:
    """
    Class to hold logical dependencies for docker images:
      git_dependencies: glob or iterable of globs that specify paths in the repo (relative to project folder). If any
                        file matching the glob(s) changes, that indicates the image needs to be rebuilt.
      docker_dependencies: mapping from docker image to input argument. The keys of the mapping indicate docker images
                           that this dockerfile derives from, and the values indicate corresponding ARG values that
                           must be passed to docker_build to specify the correct repo/image:tag
    """
    def __init__(
            self,
            git_dependencies: Union[str, Iterable[str]],
            docker_dependencies: Mapping[str, str] = MappingProxyType({})
    ):
        self.git_dependencies = git_dependencies
        self.docker_dependencies = docker_dependencies

    def has_change(self, changed_files: List[str]) -> bool:
        if isinstance(self.git_dependencies, str):
            return any(fnmatch.filter(changed_files, self.git_dependencies))
        else:
            return any(match
                       for dependency in self.git_dependencies
                       for match in fnmatch.filter(changed_files, dependency))

    def depends_on(self, changed_docker_images: Set[str]) -> bool:
        return not changed_docker_images.isdisjoint(self.docker_dependencies.keys())


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
        "delly": ImageDependencies("dockerfiles/delly/*"),
        "manta": ImageDependencies("dockerfiles/manta/*"),
        "melt": ImageDependencies("dockerfiles/melt/*"),
        "wham": ImageDependencies("dockerfiles/wham/*"),
        "sv-base-mini": ImageDependencies("dockerfiles/sv-base-mini/*"),
        "samtools-cloud": ImageDependencies("dockerfiles/samtools-cloud/*", {"sv-base-mini": "MINIBASE_IMAGE"}),
        "sv-base-virtual-env": ImageDependencies("dockerfiles/sv-base-virtual-env/*"),
        "sv-base": ImageDependencies(
            "dockerfiles/sv-base/*",
            {"sv-base-mini": "MINIBASE_IMAGE", "sv-base-virtual-env": "VIRTUAL_ENV_IMAGE"}),
        "cnmops-virtual-env": ImageDependencies(
            "dockerfiles/cnmops-virtual-env/*",
            {"sv-base-virtual-env": "VIRTUAL_ENV_IMAGE"}),
        "cnmops": ImageDependencies(
            ("dockerfiles/cnmops/*", "src/WGD/*"),
            {"sv-base": "SVBASE_IMAGE", "cnmops-virtual-env": "VIRTUAL_ENV_IMAGE"}),
        "sv-pipeline-base-virtual-env": ImageDependencies(
            "dockerfiles/sv-pipeline-base-virtual-env/*",
            {"sv-base-virtual-env": "VIRTUAL_ENV_IMAGE"}),
        "sv-pipeline-base": ImageDependencies(
            ("dockerfiles/sv-pipeline-base/*", "src/*"),
            {"sv-base": "SVBASE_IMAGE", "sv-pipeline-base-virtual-env": "VIRTUAL_ENV_IMAGE"}),
        "sv-pipeline-virtual-env": ImageDependencies(
            "dockerfiles/sv-pipeline-virtual-env/*",
            {"sv-pipeline-base-virtual-env": "VIRTUAL_ENV_IMAGE"}),
        "sv-pipeline": ImageDependencies(
            "dockerfiles/sv-pipeline/*",
            {"sv-pipeline-base": "SV_PIPELINE_BASE_IMAGE", "sv-pipeline-virtual-env": "VIRTUAL_ENV_IMAGE"}),
        "sv-pipeline-hail-virtual-env": ImageDependencies(
            "dockerfiles/sv-pipeline-hail-virtual-env/*",
            {"sv-pipeline-virtual-env": "VIRTUAL_ENV_IMAGE"}),
        "sv-pipeline-hail": ImageDependencies(
            "dockerfiles/sv-pipeline-hail/*",
            {"sv-pipeline": "SV_PIPELINE_IMAGE", "sv-pipeline-hail-virtual-env": "VIRTUAL_ENV_IMAGE"}),
        "sv-pipeline-updates-virtual-env": ImageDependencies(
            "dockerfiles/sv-pipeline-updates-virtual-env/*",
            {"sv-pipeline-virtual-env": "VIRTUAL_ENV_IMAGE"}),
        "sv-pipeline-updates": ImageDependencies(
            "dockerfiles/sv-pipeline-updates/*",
            {"sv-pipeline": "SV_PIPELINE_IMAGE", "sv-pipeline-updates-virtual-env": "VIRTUAL_ENV_IMAGE"}),
        "sv-pipeline-children-virtual-env": ImageDependencies(
            "dockerfiles/sv-pipeline-children-virtual-env/*",
            {"sv-pipeline-virtual-env": "VIRTUAL_ENV_IMAGE"}),
        "sv-pipeline-qc-virtual-env": ImageDependencies(
            "dockerfiles/sv-pipeline-qc-virtual-env/*",
            {"sv-pipeline-children-virtual-env": "VIRTUAL_ENV_IMAGE"}),
        "sv-pipeline-qc": ImageDependencies(
            "dockerfiles/sv-pipeline-qc/*",
            {"sv-pipeline": "SV_PIPELINE_IMAGE", "sv-pipeline-qc-virtual-env": "VIRTUAL_ENV_IMAGE"}),
        "sv-pipeline-rdtest-virtual-env": ImageDependencies(
            "dockerfiles/sv-pipeline-rdtest-virtual-env/*",
            {"sv-pipeline-children-virtual-env": "VIRTUAL_ENV_IMAGE"}),
        "sv-pipeline-rdtest": ImageDependencies(
            "dockerfiles/sv-pipeline-rdtest/*",
            {"sv-pipeline": "SV_PIPELINE_IMAGE", "sv-pipeline-rdtest-virtual-env": "VIRTUAL_ENV_IMAGE"}),
    }
    non_public_images = frozenset({'melt'})
    images_built_by_all = frozenset(dependencies.keys()).difference({"melt"})
    accepted_target_values = frozenset(dependencies.keys()).union({"all"})
    latest_tag = "latest"

    def __init__(
            self,
            project_arguments: argparse.Namespace,
            launch_script_path: str = Paths.this_script_folder
    ):
        if project_arguments.docker_repo is not None:
            os.system("docker login")

        self.project_arguments = project_arguments
        self.launch_script_path = launch_script_path
        self.working_dir = None
        self.build_priority = {}
        self.dockers_json = ProjectBuilder.load_json(self.project_arguments.input_json)
        self.current_docker_images = {}

    @staticmethod
    def load_json(json_file: str) -> Dict[str, str]:
        if not os.path.isfile(json_file):
            return {}
        with open(json_file, 'r') as f_in:
            return json.load(f_in)

    @staticmethod
    def get_target_from_image(docker_image: str) -> str:
        # note: will process any string, even if it does not have '/' or ':' in it, even if it is not a docker image
        # (in which case it will return the whole docker_image string)
        return docker_image.rsplit('/', 1)[-1].split(':', 1)[0]

    @staticmethod
    def get_image_repo(docker_image: str) -> Optional[str]:
        repo_ind = docker_image.rfind('/')
        return None if repo_ind < 0 else docker_image[:repo_ind]

    @staticmethod
    def is_image_local(docker_image: str) -> bool:
        return ProjectBuilder.get_image_repo(docker_image) is None

    def get_current_image(self, target: str, throw_error_on_no_image: bool = True) -> str:
        current_image = self.current_docker_images.get(target, None)
        if current_image is None:
            for image in self.dockers_json.values():
                if ProjectBuilder.get_target_from_image(image) == target:
                    self.current_docker_images[target] = image
                    current_image = image
                    break
        if throw_error_on_no_image and current_image is None:
            raise ValueError(
                f"{target} has no current image. ProjectBuilder should build images in order. This is a bug!"
            )
        return current_image

    def update_dockers_json(self):
        output_json = self.project_arguments.output_json
        if output_json is None or not output_json or output_json == Paths.dev_null:
            return  # no update is desired
        new_dockers_json = {
            json_key: (self.get_current_image(ProjectBuilder.get_target_from_image(docker_image))
                       if ProjectBuilder.get_target_from_image(docker_image) in self.dependencies
                       else docker_image)

            for json_key, docker_image in self.dockers_json.items()
        }
        # if a new image has been added that is not used by dockers json, store it as a distinct value to have a record
        # of what tag is current
        dockers_json_images = set(new_dockers_json.values())
        for target, image in self.current_docker_images.items():
            if image not in dockers_json_images:
                if target in new_dockers_json:
                    # this requires coincidences bordering on malicious, but check and throw a sensible error message
                    raise ValueError(
                        f"Unable to update {output_json} because {image} is not used in input dockers json but its"
                        f"target name {target} conflicts with an existing key."
                    )
                new_dockers_json[target] = image

        old_dockers_json = ProjectBuilder.load_json(output_json)

        if new_dockers_json != old_dockers_json:
            # update dockers.json with the new data
            if self.project_arguments.dry_run:
                print(f"Write output dockers json at {output_json}")
                print(json.dumps(new_dockers_json, indent="  "))
            else:
                with open(output_json, 'w') as f_out:
                    json.dump(new_dockers_json, f_out, indent="  ")

    def get_build_priority(self, target_name: str) -> int:
        if target_name not in self.build_priority:
            build_deps = ProjectBuilder.dependencies[target_name].docker_dependencies
            self.build_priority[target_name] = 0 if not build_deps \
                else 1 + max((self.get_build_priority(build_dep) for build_dep in build_deps.keys()), default=0)
        return self.build_priority[target_name]

    def _add_image_prereqs(self, build_targets: Set[str]) -> Set[str]:
        # Ensure that image prerequisite that a target requires exists. If not, add them to targets to build.
        # (This should almost never happen unless someone has messed with dockers.json or added a brand-new image.)
        prereqs = {
            prereq
            for target in build_targets
            for prereq in self.dependencies[target].docker_dependencies.keys()
            if self.get_current_image(prereq, throw_error_on_no_image=False) is None
        }.difference(build_targets)
        while prereqs:
            build_targets.update(prereqs)
            prereqs = {
                prereq
                for target in prereqs
                for prereq in self.dependencies[target].docker_dependencies.keys()
                if self.get_current_image(prereq, throw_error_on_no_image=False) is None
            }.difference(prereqs)
        return build_targets

    def get_ordered_build_chain_list(self) -> List[str]:
        # seed with all targets that should be built
        new_targets_to_build = self._add_image_prereqs(
            {
                target for target in self.project_arguments.targets
                if not ImageBuilder(target, self).do_not_rebuild
            }
        )

        if self.project_arguments.skip_dependent_images:
            # only build targets
            targets_to_build = new_targets_to_build
        else:
            # need to build all the targets, and everything that depends on them
            # iteratively add dependent images until nothing needs to be added
            targets_to_build = set()
            while new_targets_to_build:
                targets_to_build.update(new_targets_to_build)
                new_targets_to_build = self._add_image_prereqs(
                    {
                        image
                        for image, dependencies in self.dependencies.items()
                        if dependencies.depends_on(new_targets_to_build) and
                        not ImageBuilder(image, self).do_not_rebuild
                    }.difference(targets_to_build)
                )

        # noinspection PyTypeChecker
        return sorted(targets_to_build, key=self.get_build_priority)

    def go_to_workdir(self) -> str:
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
    def remote_docker_repos(self) -> Tuple[str, ...]:
        return () if self.project_arguments.docker_repo is None \
            else (self.project_arguments.docker_repo,) if isinstance(self.project_arguments.docker_repo, str) \
            else tuple(self.project_arguments.docker_repo)

    @property
    def built_images(self) -> List[str]:
        # note: need to separate out format part of string because docker format conflicts with f-string
        return get_command_output(
            f"docker images *:{self.project_arguments.image_tag} " + "--format '{{.Repository}}:{{.Tag}}'"
        ).split()

    def build_and_push(self):
        print("Project args:")
        pprint.pprint(vars(self.project_arguments))
        print("")
        possible_tmp_dir_path = self.go_to_workdir()
        build_completed = False

        if os.system("docker system info > /dev/null 2>&1") != 0:
            raise RuntimeError("docker daemon is not running or returned error")

        try:
            # set the appropriate targets
            if self.project_arguments.base_git_commit is not None:
                changed_project_files = self.changed_project_files
                print(f"changed_project_files: {changed_project_files}", file=sys.stderr)
                self.project_arguments.targets = [
                    target for target, image_dependencies in self.dependencies.items()
                    if image_dependencies.has_change(changed_project_files)
                ]
                print(f"targets = {self.project_arguments.targets}")
            elif "all" in self.project_arguments.targets:
                self.project_arguments.targets = ProjectBuilder.images_built_by_all
            # get targets + their dependencies in order (so that builds succeed)
            expanded_build_targets = self.get_ordered_build_chain_list()

            # build each required dependency
            print("Building and pushing the following targets in order:" if self.remote_docker_repos else
                  "Building the following targets in order:")
            print(expanded_build_targets)
            if self.project_arguments.dry_run:
                for target_name in expanded_build_targets:
                    ImageBuilder(target_name, self).update_current_image()
            else:
                print(colored('#' * 50, 'magenta'))
                for target_name in expanded_build_targets:
                    a = colored("Building image ", "grey")
                    b = colored(target_name + ":" + self.project_arguments.image_tag, "yellow", attrs=['bold'])
                    c = colored(" ...", "grey")
                    print(a, b, c)

                    build_time_args = {
                        arg: self.get_current_image(image_name)
                        for image_name, arg in ProjectBuilder.dependencies[target_name].docker_dependencies.items()
                    }

                    image_builder = ImageBuilder(target_name, self)
                    image_builder.build(build_time_args)
                    image_builder.push()
                    print(colored('#' * 50, 'magenta'))

                print(colored('BUILD PROCESS SUCCESS!', 'green'))

            if self.remote_docker_repos:
                # push any images that are purely local (this can happen if e.g. images are built without specifying
                # --docker-repo during development, but upon successful build the developer wants to push them)
                local_images_to_push = [
                    image
                    for image in [
                        self.get_current_image(target, throw_error_on_no_image=False)
                        for target in self.dependencies.keys()
                    ]
                    if image is not None and ProjectBuilder.is_image_local(image)
                ]
                if local_images_to_push:
                    print(f"Also push local images {local_images_to_push} to {self.remote_docker_repos}")
                if not self.project_arguments.dry_run:
                    for local_image in local_images_to_push:
                        ImageBuilder(local_image, self).push()

            build_completed = True
        finally:
            self.update_dockers_json()
            if not self.project_arguments.skip_cleanup:
                self.cleanup(possible_tmp_dir_path, build_completed)

    def cleanup(self, tmp_dir_path: str, build_completed):
        if not self.project_arguments.dry_run:  # don't clean up docker stuff on dry-run
            if build_completed:
                # clean dangling images (i.e. those "<none>" images), stopped containers, etc
                os.system("docker system prune -f")
            else:
                # clean stopped containers, but leave images so that cached layers go faster next time
                os.system("docker container prune -f")

        # "rm -rf" staging dir, if was specified
        if (self.project_arguments.staging_dir is not None) and (tmp_dir_path is not None):
            os.chdir(self.launch_script_path)  # first cd back to launch_script_path
            shutil.rmtree(tmp_dir_path)

    @property
    def changed_project_files(self) -> List[str]:
        base_git_commit = self.project_arguments.base_git_commit
        current_git_commit = self.project_arguments.current_git_commit
        return get_command_output(
            f"git diff --name-only {base_git_commit}" if current_git_commit is None else
            f"git diff --name-only {base_git_commit} {current_git_commit}"
        ).strip().split('\n')


class ImageBuilder:  # class for building and pushing a single image
    def __init__(self, name: str, project_builder: ProjectBuilder):
        if ':' in name:
            self.name, self.tag = name.split(':', 1)
        else:
            self.name = name
            self.tag = project_builder.project_arguments.image_tag
        self.project_builder = project_builder

    @property
    def working_dir(self) -> str:
        return self.project_builder.working_dir

    @property
    def local_image(self) -> str:
        return f"{self.name}:{self.tag}"

    @property
    def remote_docker_repos(self) -> Tuple[str, ...]:
        return self.project_builder.remote_docker_repos

    @property
    def remote_images(self) -> Tuple[str, ...]:
        remote_tags = (self.tag, ProjectBuilder.latest_tag) if self.project_builder.project_arguments.update_latest \
            else (self.tag,)
        return tuple(
            f"{repo}/{self.name}:{remote_tag}"
            for repo in self.remote_docker_repos
            for remote_tag in remote_tags
        )

    @property
    def no_force_rebuild(self) -> bool:
        return self.project_builder.project_arguments.no_force_rebuild

    @staticmethod
    def remote_image_exists(remote_image: str) -> bool:
        output, stderr, return_code = get_command_output(
            f"DOCKER_CLI_EXPERIMENTAL=enabled docker manifest inspect {remote_image}",
            raise_on_error=False,
            return_error_info=True
        )
        if return_code != 0 and "Caller does not have permission" in stderr:
            # non-public images (e.g. melt) may be stored in restrictive repos. Return True because we don't want to
            # mess with them
            return True
        else:
            return return_code == 0

    @staticmethod
    def image_is_built(local_image: str) -> bool:
        images = {
            image for image in get_command_output(
                'docker images --format "{{.Repository}}:{{.Tag}}"'
            ).split('\n')
        }
        return local_image in images

    @property
    def do_not_rebuild(self) -> bool:
        return self.no_force_rebuild and ImageBuilder.image_is_built(self.local_image)

    def build(self, build_time_args: Mapping[str, str]):
        if self.do_not_rebuild:
            print(f"skipping build of {self.local_image} because it is already built and --no-force-rebuild is set")
            return
        # standard build command
        docker_build_command = "docker build --progress plain \\\n    "
        docker_build_command += "-f " + f"{self.working_dir}/dockerfiles/{self.name}/Dockerfile" + " \\\n    "
        docker_build_command += "--tag " + self.local_image + " \\\n    "
        # parse extra args list
        for key, value in build_time_args.items():
            docker_build_command += "--build-arg " + key + "=" + value + " \\\n    "

        will_push = 0 != len(self.remote_docker_repos) and any(e is not None for e in self.remote_docker_repos)
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
        """
        Push everything related to this image to remote repo
        """
        for remote_image in self.remote_images:
            # do not push images with very restrictive licenses
            if self.name in ProjectBuilder.non_public_images and not remote_image.startswith('us.gcr.io'):
                print(colored(f"Refusing to push non-public image {self.name} to non us.grc.io repo at {remote_image}",
                              "red"))
                continue
            if ImageBuilder.docker_tag(self.local_image, remote_image) != 0:
                raise RuntimeError(f"Failed to tag image ({remote_image}) for pushing to remote")
            if ImageBuilder.docker_push(remote_image) != 0:
                raise RuntimeError(f"Failed to push image {remote_image}")

        self.update_current_image()

    def update_current_image(self):
        self.project_builder.current_docker_images[self.name] = \
            self.remote_images[0] if self.remote_docker_repos else self.local_image

    @staticmethod
    def docker_tag(source_image: str, target_image: str) -> int:
        tag_command = f"docker tag {source_image} {target_image}"
        print(tag_command)
        return os.system(tag_command)

    @staticmethod
    def docker_push(remote_image: str) -> int:
        docker_push_command = f"docker push {remote_image}"
        print(docker_push_command)
        return os.system(docker_push_command)


def get_command_output(
        command: str,
        encoding: str = "utf-8",
        raise_on_error: bool = True,
        return_error_info: bool = False
) -> Union[str, Tuple[str, str, int]]:
    """
    Execute shell command. Raise exception if unsuccessful, otherwise return string with output
    """
    sub_p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    with sub_p.stdout as pipe_in, sub_p.stderr as pipe_err:
        output = pipe_in.read().decode(encoding)
        stderr = pipe_err.read().decode(encoding)
    return_code = sub_p.poll()
    if raise_on_error and return_code != 0:
        raise RuntimeError('Error executing %s:\n%s' % (command, stderr[:-1]))
    return (output, stderr, return_code) if return_error_info else output


class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    pass


def __parse_arguments(args_list: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Tool for building docker images for GATK-SV pipeline" + "\n" + program_operation_description,
        formatter_class=CustomFormatter,
    )
    # required arguments
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
    docker_remote_args_group.add_argument('--docker-repo', type=str,
                                          help='Docker repo to push images to. This will push images that are built '
                                               'this run of build_docker.py, or that currently have only a local image '
                                               'in --input-json')
    docker_remote_args_group.add_argument('--gcr-project', type=str,
                                          help='Deprecated. Used to determine which docker repo to push images to. Use '
                                               '--docker-repo instead.')
    docker_remote_args_group.add_argument('--update-latest', action='store_true',
                                          help=f'also update \"{ProjectBuilder.latest_tag}\" tag in remote docker'
                                               f'repo(s)')
    docker_remote_args_group.add_argument('--input-json', type=str, default=Paths.dockers_json_path,
                                          help="Path to dockers.json to use as input. This file serves as a store for "
                                               "both the default docker image to use for various gatk-sv WDLs, and for "
                                               "the most up-to-date docker tag for each docker image.")
    docker_remote_args_group.add_argument('--output-json', type=str, default=Paths.dockers_json_path,
                                          help=f"Path to output updated dockers.json. Set to {Paths.dev_null} to turn "
                                               "off updates")
    parser.add_argument("--dry-run", action="store_true",
                        help="Compute docker images that will be build, but don't actually build or push")
    parser.add_argument("--skip-dependent-images", action="store_true",
                        help="Don't build images that depend on targets. Can leave images in an unreproducible state: "
                             "only use this if you are sure you know what you are doing!")
    parser.add_argument(
        '--targets', nargs='*', type=str,
        help='Manually-specified list project docker image(s) you want to build (note "all" does not include melt). '
             'Alternatively can specify --base-git-commit/--current-git-commit to automatically determine targets.'
    )
    short_git_hash_head = get_command_output("git rev-parse --short HEAD").strip()
    parser.add_argument('--image-tag', type=str, default=short_git_hash_head,
                        help='tag to be applied to all images being built')
    parser.add_argument(
        # flag to turn off git protection (default mode is refusing to build when there are untracked files and/or
        # uncommitted changes)
        '--disable-git-protect', action='store_true',
        help='disable git check/protect when building from local files (will use uncommited changes to build)'
    )
    parser.add_argument('--no-force-rebuild', action='store_true',
                        help='Do not rebuild docker images if the exact image and tag already exist.')
    parser.add_argument('--skip-cleanup', action='store_true',
                        help='skip cleanup after successful and unsuccessful build attempts.')
    parser.add_argument('--base-git-commit', type=str,
                        help="This script can have targets specified manually (via --targets) or it can automatically "
                             "determine which docker image \"targets\" to build by examining which files have changed "
                             "in the git repo. When auto-determining build targets, this options specifies the baseline"
                             " git commit to check for changes, i.e. only files altered since this commit should be "
                             "considered changed. Can be a SHA or other  specifier "
                             "(e.g. HEAD^)")
    parser.add_argument('--current-git-commit', type=str,
                        help="This script can automatically determine which docker image \"targets\" to build by "
                             "examining which files have changed in the git repo. When auto-determining build targets, "
                             "this options specifies the current git commit to check for changes. If omitted,"
                             " use current status of git repo with uncommitted changes. Can be a SHA or other specifier"
                             " (e.g. HEAD)")

    # parse and consistency check
    if len(args_list) <= 1:  # no arguments, print help and exit with success
        parser.parse_args(['-h'])
        sys.exit(0)
    parsed_args = parser.parse_args(args_list[1:])

    if not os.access(parsed_args.input_json, os.R_OK):
        raise ValueError("--input-json must specify a path to a file with read access.")

    if parsed_args.base_git_commit is None:
        if parsed_args.targets is None:
            raise ValueError('Must specify exactly one of "--base-git-commit" or "--targets", but neither were passed')
        # if passed targets are in the accepted values
        for tar in parsed_args.targets:
            if tar not in ProjectBuilder.accepted_target_values:
                raise ValueError("\"" + tar + "\" not in allowed target values")

        if "all" in parsed_args.targets:
            if 1 != len(parsed_args.targets):
                raise ValueError("when \"all\" is provided, no other target values allowed")
    else:
        if parsed_args.targets is not None:
            raise ValueError('Must specify exactly one of "--base-git-commit" or "--targets", but both were passed')

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
            raise ValueError(f"publishing \"{ProjectBuilder.latest_tag}\" docker images requires a remote Github "
                             "release tag")

    # if there are uncommitted changes when building from local files, raise exception
    if parsed_args.staging_dir is None and not parsed_args.disable_git_protect:
        s = os.popen("git status -s | wc -l | tr -d ' ' | tr -d '\n'").read()
        ret = int(s)
        if 0 != ret:
            raise ValueError(
                "Current directory has uncommitted changes or untracked files. Cautiously refusing to proceed."
            )

    if parsed_args.gcr_project is not None:
        if parsed_args.docker_repo is not None:
            raise ValueError("Both --gcr-project and --docker-repo were specified, but only one is allowed.")
        print(colored("--gcr-project is deprecated, use --docker-repo instead.", "red"))
        parsed_args.docker_repo = f"us.gcr.io/{parsed_args.gcr_project.strip('/')}"
    return parsed_args


def main(arguments: List[str]):
    project_args = __parse_arguments(arguments)
    ProjectBuilder(project_args).build_and_push()


if __name__ == "__main__":
    main(sys.argv)
