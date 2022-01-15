"""
The intent of this script is to update `dockers.json` with
the images rebuilt and pushed by `build_docker.py`. Since for
a given list of images, `build_docker.py` may rebuild additional
images, the list of all the built and pushed images may be different
from the user-requested list. However, the `build_docker.py` does
not export the list of all the built images, hence, this script
takes a bruteforce approach and checks for all the images if they
are updated by a given tag (i.e., the tag used by the `build_docker.py`
to tag all the updated images).

This script is mainly developed to be used in Github actions
for the CI/CD pipeline, and is a temporary solution until
`build_docker.py` is updated to export the list of all the built and
pushed images.
"""

import argparse
import copy
import json
import os
from enum import Enum
from os.path import dirname, abspath, join
from subprocess import Popen, PIPE
from typing import Dict, Optional, List


COLOR_RED = "\033[91m"
COLOR_GREEN = "\033[92m"
COLOR_YELLOW = "\033[93m"
COLOR_ENDC = "\033[0m"


class Status(Enum):
    NO_MARK = 0
    NEUTRAL = 1
    SUCCESS = 2
    FAILED = 3
    SKIPPED = 4


def pprint_log(msg: str, status=Status.NO_MARK, line_break=True):
    if status != Status.NO_MARK:
        # "\N{<unicode name of a character>}"
        # See https://unicode.org for other characters.
        if status == Status.NEUTRAL:
            mark = "\N{black star}"
        elif status == Status.SUCCESS:
            mark = "\N{heavy check mark}"
        elif status == status.SKIPPED:
            mark = "\N{em dash}"
        else:
            mark = "\N{heavy ballot x}"
        msg = f"{mark}\t{msg}"

    if status in [Status.SUCCESS, Status.FAILED, Status.SKIPPED]:
        if status == status.SUCCESS:
            txt_color = COLOR_GREEN
        elif status == status.FAILED:
            txt_color = COLOR_RED
        else:
            txt_color = COLOR_YELLOW
        msg = f"{txt_color}{msg}{COLOR_ENDC}"
    print(msg, end="\n" if line_break else "", flush=True)


def check_image_exist(image: str) -> (bool, Optional[str]):
    """
    Checks if a given image exists; implemented based on the
    following SO answer: https://stackoverflow.com/a/52077346/947889

    This method leverages `docker` instead of a curl-based method,
    because if the method requires any authentication to query the
    docker registry, the authentication has happened in the Github
    action's context.
    """
    docker_cmd = f"docker manifest inspect {image} > /dev/null; echo $?"
    proc = Popen(docker_cmd, stdout=PIPE, stderr=PIPE, shell=True)
    (out, error) = proc.communicate()
    out = out.decode("utf-8").strip()
    error = error.decode("utf-8").strip()
    if out == "0":
        return True, None
    else:
        return False, error


def get_updated_images(
        images: Dict[str, str], tag: str, repo: str,
        exclude_images: List[str]) -> Dict[str, str]:
    # not modifying the original in case any
    # ref to the original comes in handy.
    updated_images = copy.deepcopy(images)

    # Some decorative and logging setup.
    c = 0
    # `-1` to account for `"name" : "dockers",` which is not an image.
    keys_count = len(updated_images) - 1

    pprint_log(f"Checking {keys_count} images with the tag `{tag}`.",
               line_break=True)
    pprint_log("The result of asserting every image in the provided JSON file "
               "will be reported according to the following legend.\n",
               line_break=True)
    pprint_log("The image tag is identical to the tag in the provided JSON "
               "file, OR the tag is different but an image with the given "
               "tag does NOT exist in the container registry; hence "
               "the image listed in the JSON file will remain unchanged.\n",
               status=Status.NEUTRAL, line_break=True)
    pprint_log("The image tag is different from the tag in the provided "
               "JSON file, and an image with the given tag exists in the "
               "container registry; hence the provided JSON file will be "
               "updated to reference this image.\n",
               status=Status.SUCCESS, line_break=True)
    pprint_log("There was an error querying the image "
               "from the image registry.\n",
               status=Status.FAILED, line_break=True)
    pprint_log("The image is skipped as specified by the "
               "`--exclude-images` argument.\n",
               status=Status.SKIPPED, line_break=True)

    for name, image in updated_images.items():
        if name == "name":
            continue

        c += 1
        pprint_log(f"[{c}/{keys_count}]\t", Status.NO_MARK, line_break=False)

        if name in exclude_images:
            pprint_log(name, Status.SKIPPED)
            continue

        image_base = image.split(":")[0]
        expected_image = f"{image_base}:{tag}"
        if expected_image == image:
            pprint_log(name, Status.NEUTRAL)
            continue

        def update_existing_image():
            updated_images[name] = expected_image
            pprint_log(name, Status.SUCCESS)

        exists, error = check_image_exist(expected_image)
        if exists:
            update_existing_image()
        else:
            image_name = image_base.split("/")[-1]
            expected_image = f"{repo}/{image_name}:{tag}"
            exists, error = check_image_exist(expected_image)
            if exists:
                update_existing_image()
            else:
                # There might be authentication issues checking if the image
                # exist; hence checking for the original image's existence;
                # if this fails too, then there are issues unrelated to the
                # given tag, else we can assume an image with the given tag
                # does not exist (in other words, the image is not updated).
                e, _ = check_image_exist(image)
                if not e:
                    pprint_log(f"{name}\t{error}", Status.FAILED)
                else:
                    pprint_log(name, Status.NEUTRAL)
                    continue
    return updated_images


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="For a given Docker image tag, this script checks "
                    "if any of the Docker images listed in the provided JSON "
                    "file (e.g., `input_values/dockers.json`) are pushed to "
                    "the container registry (e.g., Google Container "
                    "Registry) with that tag.\n\n"
                    "If the given tag is different from the tag in the "
                    "provided JSON file, and an image with the given tag "
                    "exists in the registry, the script updates the "
                    "referenced Docker image; otherwise, it leaves the "
                    "referenced image unchanged. The script persists the "
                    "updated list of images in the provided output JSON file."
                    "\n\n"
                    "This script assumes an image with updated tag is "
                    "stored under the same registry as the previous one. "
                    "For example, \"us.gcr.io/repo/image:abc\" is expected "
                    "to update to \"us.gcr.io/repo/image:def\" where the "
                    "registry (i.e., `us.gcr.io`) and repository "
                    "(i.e., `repo`) are not changed."
                    "\n\n"
                    "The intended use-case of this script is for a scenario "
                    "when you know the tag used to build and push GATK-SV "
                    "Docker images using the `build_docker.py` script, "
                    "however, you're not sure which images are updated "
                    "(since `build_docker.py` may rebuild and push additional "
                    "images to the specified images, e.g. its dependencies). "
                    "This script is predominantly developed for the Github "
                    "actions, and is a temporary solution until the "
                    "`build_docker.py` is re-written to provide this "
                    "functionality.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        "image_tag",
        help="Docker image tag.")

    parser.add_argument(
        "docker_repo",
        help="Sets the repository expected to contain "
             "the images with the given tag. The repository should be fully "
             "quantifying; e.g., `us.gcr.io/broad-dsde-methods/gatk-sv`."
    )

    parser.add_argument(
        "-i", "--input-json",
        default=join(dirname(dirname(dirname(abspath(__file__)))),
                     "input_values/dockers.json"),
        help="A JSON file containing Docker images to check.")

    parser.add_argument(
        "-o", "--output-json",
        help="A JSON file to persist a list of the updated Docker images. "
             "Defaults to the input `--input-json`. Note that if the given "
             "file exists, this script will replace it."
    )

    parser.add_argument(
        "-e", "--exclude-images",
        required=False,
        nargs='+',
        default=["gatk_docker", "gatk_docker_pesr_override",
                 "genomes_in_the_cloud_docker", "linux_docker",
                 "cloud_sdk_docker"],
        help="Sets a list of docker images to be ignored. Defaults to the "
             "list of docker images not currently built by `build_docker.py`.")

    return parser.parse_args()


def main():
    args = parse_arguments()

    # Checks if the file exists, if not, errors-out.
    if not os.path.isfile(args.input_json):
        raise ValueError(f"dockers_json {args.input_json} is not a file")

    # Loads the JSON object from `args.dockers_json`.
    # Errors-out if an invalid JSON is provide.
    with open(args.input_json, "r") as f:
        images = json.load(f)

    ouput_json = args.output_json \
        if args.output_json else args.input_json

    if not os.access(os.path.dirname(ouput_json), os.W_OK):
        raise OSError(f"Unable to write to updated dockers folder "
                      f"{os.path.dirname(ouput_json)}")

    updated_images = get_updated_images(
        images, args.image_tag, args.docker_repo.strip("/"),
        args.exclude_images)
    with open(ouput_json, "w") as f:
        json.dump(updated_images, f, indent=2)


if __name__ == '__main__':
    main()
