import argparse, os, os.path
import time
import tempfile, shutil
from termcolor import colored

############################################################

# WARN: do not re-order this trivially, because
#       except the first row off-the-shelf callers, the rest of dockers have non-trivial dependency structures
INDIVIDUAL_TARGET_VALUES = ("delly", "manta", "melt", "wham",
                            "sv-base-mini",
                            "sv-base", "samtools-cloud",
                            "cnmops", "sv-pipeline-base", "sv-pipeline-children-r",
                            "sv-pipeline", "sv-pipeline-rdtest", "sv-pipeline-qc")
ACCEPTED_TARGET_VALUES = ("all",) + INDIVIDUAL_TARGET_VALUES

############################################################
# dummy Exceptions for use when OS docker build fails
class UserError(Exception):
    pass

class DockerBuildError(Exception):
    pass

############################################################
# parsing and checking arguments
class CMD_line_args_parser:

    def __init__(self, args_list):
        parser = argparse.ArgumentParser(description='='*50 + "\nBuilding docker images for GATK-SV pipeline v1.\n" + '='*50,
                                         formatter_class=argparse.RawDescriptionHelpFormatter)

        # required arguments
        required_args_group = parser.add_argument_group('Required', 'required arguments')

        required_args_group.add_argument('--targets',
                                         nargs    = '+',
                                         type     = str,
                                         required = True,
                                         help     = 'the sub project docker(s) you want to build (note "all" does not include melt)')

        required_args_group.add_argument('--image-tag',
                                         type     = str,
                                         required = True,
                                         help     = 'tag to be applied to all images being built')

        # to build from local or remote git tag/hash values
        git_args_group = parser.add_argument_group('Mutex args', 'remote git tag/hash values (mutually exclusive)')
        git_mutex_args_group = git_args_group.add_mutually_exclusive_group()

        git_mutex_args_group.add_argument('--remote-git-tag',
                                          type = str,
                                          help = 'release tag on Github; this indicates pulling from Github to a staging dir')

        git_mutex_args_group.add_argument('--remote-git-hash',
                                          type = str,
                                          help = 'a hash value on Github; this indicates pulling from Github to a staging dir')

        # to build from remote/Github (staging required)
        remote_git_args_group = parser.add_argument_group('Remote git',
                                                          'args involved when building from remote git tags/hashes')

        remote_git_args_group.add_argument('--staging-dir',
                                           type = str,
                                           help = 'a temporary staging directory to store builds; required only when pulling from Github, ignored otherwise')

        remote_git_args_group.add_argument('--use-ssh',
                                           action = 'store_true',
                                           help   = 'use SSH to pull from github')

        # flag to turn on push to Dockerhub and/or GCR
        docker_remote_args_group = parser.add_argument_group('Docker push',
                                                             'controlling behavior related pushing dockers to remote repos')

        docker_remote_args_group.add_argument('--gcr-project',
                                              type = str,
                                              help = 'GCR billing project to push the images to')

        docker_remote_args_group.add_argument('--update-latest',
                                              action = 'store_true',
                                              help   = 'also update \"latest\" tag in remote docker repo(s)')

        # flag to turn off git protection (default modue is refusing to build when there are untracked files and/or uncommited changes)
        parser.add_argument('--disable-git-protect',
                            action = 'store_true',
                            help = 'disable git check/protect when building from local files (will use uncommited changes to build)')

        parser.add_argument('--skip-base-image-build',
                            action = 'store_true',
                            help = 'skip rebuild of the target\'s base image(s). Assumes that the base image(s) already exist with same tag.')

        parser.add_argument('--skip-cleanup',
                            action = 'store_true',
                            help = 'skip cleanup after successful and unsuccessful build attempts. This will speed up subsequent builds.')

        # parse and consistency check
        parsed_args = parser.parse_args(args_list)
        CMD_line_args_parser.consistency_check(parsed_args)
        self.project_args = parsed_args


    # cmd line args consistency check
    @staticmethod
    def consistency_check(argparse_namespace_obj):

        # if passed targets are in the accepted values
        for tar in argparse_namespace_obj.targets:
            if(tar not in ACCEPTED_TARGET_VALUES):
                raise UserError("\"" + tar + "\" not in allowed target values")

        if ("all" in argparse_namespace_obj.targets):
            if 1 != len(argparse_namespace_obj.targets):
                raise UserError("when \"all\" is provided, no other target values allowed")

        # if "use_ssh" flag is turned on, remote git tag/hash should be provided
        if (argparse_namespace_obj.use_ssh is True):
            if (argparse_namespace_obj.remote_git_tag is None) and (argparse_namespace_obj.remote_git_hash is None):
                raise UserError("\"use_ssh\" is specified but remote git tag/hash is not")

        # if remote git tag/hash and/or is specified, staging dir should be specified
        if (argparse_namespace_obj.remote_git_tag is not None) or (argparse_namespace_obj.remote_git_hash is not None):
            if argparse_namespace_obj.staging_dir is None:
                raise UserError("remote git tag/hash is specified but staging_dir is not")

        # if requesting to update "latest" tag in remote docker repo(s), remote git release tag must be specified
        if (argparse_namespace_obj.update_latest is True):
            if argparse_namespace_obj.remote_git_tag is None:
                raise UserError("publishing \"latest\" docker images requires a remote Github release tag")

        # if there're un-committed changes when building from local files, raise exception
        if (argparse_namespace_obj.staging_dir is None and not argparse_namespace_obj.disable_git_protect):
            s = os.popen("git status -s | wc -l | tr -d ' ' | tr -d '\n'").read()
            ret = int(s)
            if 0 != ret:
                raise UserError("Current directory has uncommited changes or untracked files. Cautiously refusing to proceed.")

############################################################
# controlling the build and push of a single image
class Docker_Build:

    NON_PUBLIC_DOCKERS = ('melt')

    def __init__(self, name, tag, build_context, remote_docker_repos):

        self.name = name
        self.tag = tag
        self.build_context = build_context
        self.remote_docker_repos = remote_docker_repos

    def build(self, built_time_args_dict):

        # get to the requested directory
        docker_build_command  = "cd " + self.build_context + " && \\\n"
        # standard build command
        docker_build_command += "docker build --progress plain \\\n    "
        docker_build_command += "--tag " + self.name + ":" + self.tag + " \\\n    "
        # parse extra args list
        for key, value in built_time_args_dict.items():
            docker_build_command += "--build-arg " + key + "=" + value + " \\\n    "

        will_push = 0!=len(self.remote_docker_repos) and any(e is not None for e in self.remote_docker_repos)
        docker_build_command += "--squash . " if (will_push) else ". "

        # build and time it
        print(docker_build_command)
        start_time = time.time()
        ret = os.system(docker_build_command)
        if 0 != ret:
            raise DockerBuildError("Failed to build image " + self.name + ":" + self.tag)
        elapsed_time = time.time() - start_time
        elapsed_min, elapsed_sec = divmod(elapsed_time, 60)
        print("Time spent on docker build:")
        print(f"{elapsed_min} minutes, {elapsed_sec} seconds")

    def push(self, is_update_latest):

        for rep in self.remote_docker_repos:
            # do not push images with very restrictive licenses
            if (self.name in Docker_Build.NON_PUBLIC_DOCKERS) and ( not rep.startswith('us.gcr.io') ):
                print(colored("Refusing to push non-public image " + self.name + " to " + rep, "red"))
                next

            remote_tag = "latest" if (is_update_latest) else self.tag
            docker_tag_command  = "docker tag " + self.name + ":" + self.tag + " " + rep + "/" + self.name + ":" + remote_tag
            docker_push_command = "docker push " + rep + "/" + self.name + ":" + remote_tag
            print(docker_tag_command)
            print(docker_push_command)
            ret = os.system(docker_tag_command)
            if 0 != ret:
                raise DockerBuildError("Failed to tag image for pushing to remote")
            ret = os.system(docker_push_command)
            if 0 != ret:
                raise DockerBuildError("Failed to push image")

############################################################
# controlling the build and push of all requested images
class Project_Build:

    INTERM_RESOURCE_IMG = 'gatksv-pipeline-v1-resources'

    GITHUB_ORG     = 'broadinstitute'
    GITHUB_REPO    = 'gatk-sv'

    #### for constructing an ordered build chain, to resolve dependency
    DEP_DICT = {'delly': None, 'manta': None, 'melt': None, 'wham': None,
                'sv-base-mini': None,
                'samtools-cloud': 'sv-base-mini',
                'sv-base': 'sv-base-mini',
                'cnmops': 'sv-base',
                'sv-pipeline-base': 'sv-base',
                'sv-pipeline': 'sv-pipeline-base',
                'sv-pipeline-children-r': 'sv-pipeline-base',
                'sv-pipeline-rdtest': 'sv-pipeline-children-r', 'sv-pipeline-qc': 'sv-pipeline-children-r'}
    BUILD_PRIORITY = {'delly': 0, 'manta': 0, 'melt': 0, 'wham': 0,
                      'sv-base-mini': 0,
                      'sv-base': 1, 'samtools-cloud': 1,
                      'cnmops': 2, 'sv-pipeline-base': 2,
                      'sv-pipeline': 3,
                      'sv-pipeline-children-r': 3,
                      'sv-pipeline-rdtest': 4, 'sv-pipeline-qc': 4}
    # for use when a single target is to be built
    @staticmethod
    def get_ordered_build_chain_single(target_name):
        chain = ()
        t = target_name
        while ( Project_Build.DEP_DICT[t] is not None):
            t = Project_Build.DEP_DICT[t]
            chain = (t,) + chain
        return chain + (target_name,)
    # for use when multiple targets are to be built
    @staticmethod
    def get_ordered_build_chain_list(target_name_list):
        dup_chain = []
        for t in target_name_list:
            dup_chain.extend( list(Project_Build.get_ordered_build_chain_single(t)) )
        agg_chain = list(set(dup_chain))  # uniquify
        build_chain = sorted(agg_chain, key=lambda name : Project_Build.BUILD_PRIORITY[name])
        return tuple(build_chain) # immutable, at least attempt



    def __init__(self, project_arguments, launch_script_path):

        # Todo: we can also check auth to push to GCR
        if project_arguments.gcr_project is not None:
            os.system("docker login")

        self.project_arguments = project_arguments
        self.launch_script_path = launch_script_path
        self.working_dir = None
        self.successfully_built_images = []


    def go_to_workdir(self):

        tmp_dir_path = None

        # local mode
        if(self.project_arguments.staging_dir is None):
            par = os.path.dirname(self.launch_script_path)
            par = os.path.dirname(par)
            self.working_dir = os.path.dirname(par)
            os.chdir(self.working_dir)

        else:
            # if staging is required, mkdir, cd, and pull
            if self.project_arguments.staging_dir.endswith("/"):
                tmp_dir_path = tempfile.mkdtemp(prefix = self.project_arguments.staging_dir)
            else:
                tmp_dir_path = tempfile.mkdtemp(prefix = self.project_arguments.staging_dir + "/")
            connect_mode = "git@github.com:" if self.project_arguments.use_ssh else "https://github.com"
            ret = os.system("git clone " +
                            connect_mode +
                            "/" + Project_Build.GITHUB_ORG + "/" + Project_Build.GITHUB_REPO + ".git " +
                            tmp_dir_path)
            if 0 != ret:
                raise UserError("Failed to clone the repo.")
            self.working_dir = tmp_dir_path
            os.chdir(self.working_dir)

        # checkout desired hash or tag, if building remotely
        if (self.project_arguments.remote_git_tag is not None):
            git_checkout_cmd = "git checkout tags/" + self.project_arguments.remote_git_tag
            ret = os.system( git_checkout_cmd )
            if 0 != ret:
                raise UserError("Seems that the provided git tag ["
                                + self.project_arguments.remote_git_tag
                                + "] does not exist")
        elif (self.project_arguments.remote_git_hash is not None):
            git_checkout_cmd = "git checkout      " + self.project_arguments.remote_git_hash
            ret = os.system( git_checkout_cmd )
            if 0 != ret:
                raise UserError("Seems that the provided git hash ["
                                + self.project_arguments.remote_git_hash
                                + "] does not exist")

        print("Working directory: " + os.getcwd())
        return tmp_dir_path

    def build_and_push(self):

        # start docker daemon, if one hasn't been started yet
        os.system("open --background -a Docker && while ! docker system info > /dev/null 2>&1; do sleep 1; done")

        # prepare resources docker
        print(colored('#################################################', 'magenta'))
        a = colored("Building intermediate resource image", "grey")
        b = colored(Project_Build.INTERM_RESOURCE_IMG + ":" + self.project_arguments.image_tag, "yellow", attrs=['bold'])
        c = colored(" ...", "grey")
        print(a, b, c)
        resource_docker_build_cmd = "docker build -t " + Project_Build.INTERM_RESOURCE_IMG + ":" + self.project_arguments.image_tag + " -f scripts/docker/resources.Dockerfile ."
        print(resource_docker_build_cmd)
        os.system(resource_docker_build_cmd)
        print(colored('#################################################', 'magenta'))

        # if build all, easy
        # otherwise construct build chain on the fly
        if "all" in self.project_arguments.targets:
            expanded_build_targets = tuple(target for target in INDIVIDUAL_TARGET_VALUES if target != "melt")
        elif self.project_arguments.skip_base_image_build:
            expanded_build_targets = self.project_arguments.targets
        else:
            expanded_build_targets = Project_Build.get_ordered_build_chain_list(tuple(self.project_arguments.targets))
        print("Building the following targets in order:")
        print(expanded_build_targets)
        print(colored('#################################################', 'magenta'))

        for proj in expanded_build_targets:

            a = colored("Building image ", "grey")
            b = colored(proj + ":" + self.project_arguments.image_tag, "yellow", attrs=['bold'])
            c = colored(" ...", "grey")
            print(a, b, c)

            build_time_args = {}
            if (proj == "sv-base" or proj == "samtools-cloud"):
                build_time_args = {
                    "MINIBASE_IMAGE" : "sv-base-mini:" + self.project_arguments.image_tag
                }
            elif (proj == "cnmops" or proj == "sv-pipeline-base"):
                build_time_args = {
                    "GATKSV_PIPELINE_V1_RESOURCES_IMAGE" : Project_Build.INTERM_RESOURCE_IMG + ":" + self.project_arguments.image_tag,
                    "SVBASE_IMAGE" : "sv-base:" + self.project_arguments.image_tag
                }
            elif (proj == "sv-pipeline-children-r"):
                build_time_args = {
                    "SV_PIPELINE_BASE_IMAGE" : "sv-pipeline-base:" + self.project_arguments.image_tag
                }
            elif (proj == "sv-pipeline-base"):
                build_time_args = {
                    "GATKSV_PIPELINE_V1_RESOURCES_IMAGE" : Project_Build.INTERM_RESOURCE_IMG + ":" + self.project_arguments.image_tag,
                    "SVBASE_IMAGE" : "sv-base:" + self.project_arguments.image_tag
                }
            elif (proj == "sv-pipeline"):
                build_time_args = {
                    "GATKSV_PIPELINE_V1_RESOURCES_IMAGE" : Project_Build.INTERM_RESOURCE_IMG + ":" + self.project_arguments.image_tag,
                    "SV_PIPELINE_BASE_IMAGE" : "sv-pipeline-base:" + self.project_arguments.image_tag
                }
            elif (proj.startswith("sv-pipeline")):
                build_time_args = {
                    "SV_PIPELINE_BASE_R_IMAGE" : "sv-pipeline-children-r:" + self.project_arguments.image_tag
                }

            build_context = self.working_dir + "/dockerfiles/" + proj

            remote_docker_repos = []
            if (self.project_arguments.gcr_project is not None):
                remote_docker_repos.append("us.gcr.io/" + self.project_arguments.gcr_project)

            docker = Docker_Build(proj, self.project_arguments.image_tag, build_context,
                                  remote_docker_repos)
            docker.build(build_time_args)
            self.successfully_built_images.append(docker)
            print(colored('#################################################', 'magenta'))

        for succ in self.successfully_built_images:
            succ.push(self.project_arguments.update_latest)

        print(colored('BUILD PROCESS SUCCESS!', 'red'))

    def cleanup(self, tmp_dir_path):

        # clean dangling images (i.e. those "<none>" images)
        clean_dangling_images = "docker images --format \"{{.ID}}\"  --filter \"dangling=true\" -q --no-trunc | xargs docker rmi"
        os.system(clean_dangling_images)

        # clean intermediate image that are not to be pushed
        os.system("docker rmi --force " + Project_Build.INTERM_RESOURCE_IMG + ":" + self.project_arguments.image_tag)

        # "rm -rf" staging dir, if was specified
        if(self.project_arguments.staging_dir is not None) and (tmp_dir_path is not None):
            os.chdir(os.path.dirname(self.launch_script_path))  # first cd back
            shutil.rmtree(tmp_dir_path)

############################################################
# static function to control the build and cleanup
def parse_and_build(parsed_project_arguments, launch_script_path):

    import pprint
    print("Project args:")
    pprint.pprint(vars(parsed_project_arguments))
    print("")

    my_build_project = Project_Build(parsed_project_arguments, launch_script_path)

    possible_tmp_dir_path = my_build_project.go_to_workdir()

    try:
        my_build_project.build_and_push()
    except UserError as a:
        raise Exception("Build Process Errored due to an assertion error!!!\n"
                        + str(a))
    except DockerBuildError as d:
        raise Exception("Build Process Errored due to a docker build error!!!\n"
                        + str(d))
    finally:
        if not my_build_project.project_arguments.skip_cleanup:
            my_build_project.cleanup(possible_tmp_dir_path)

############################################################

if __name__ == "__main__":
    import sys
    if 1 == len(sys.argv):
        parser = CMD_line_args_parser(['-h'])
    else:
        project_args = CMD_line_args_parser(sys.argv[1:]).project_args
        this_script_path = os.path.realpath('__file__')

        parse_and_build(project_args, this_script_path)
