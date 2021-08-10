"""
Requirements:
- Java
- Womtool
- carrot_cli
"""

import argparse
import subprocess
from subprocess import check_output
import time
import calendar
import json
import os
import urllib.request
from urllib.error import HTTPError

# The directories where the WDLs
# and their tests are located.
WDLS_DIR_RELATIVE = "../wdl"
WDLS_DIR = "wdl"
WDLS_TEST_DIR = "wdl_test"

WOMTOOL_PATH = "~/code/cromwell/womtool-61.jar"

DEFAULT_EVAL_INPUTS_FILENAME = "eval_input_defaults.json"
DEFAULT_TEST_INPUTS_FILENAME = "test_input_defaults.json"
EVAL_INPUTS_FILENAME = "eval_input.json"
TEST_INPUTS_FILENAME = "test_input.json"


# This is the list of the currently
# supported results type by Carrot.
# Keep the list lowercase, as carrot is
# case-sensitive and requires lower-case. See:
# https://github.com/broadinstitute/carrot/blob/1b2905e634e0887a94c63743a78b1731bd3a637c/src/custom_sql_types.rs#L158-L162
SUPPORTED_RESULTS_TYPE = ["numeric", "file", "text"]


class CarrotHelper:
    def __init__(self, working_dir=".",
                 pipelines_filename=".carrot_pipelines.json",
                 repository="https://github.com/vjalili/gatk-sv",
                 branch="ehdn_unit_test",
                 email=None):
        self.working_dir = working_dir
        self.pipelines_filename = pipelines_filename
        self.repository = repository
        self.branch = branch
        self.email = email or self._get_email_from_config()
        self.pipelines = {}
        self.load_pipelines()
        if not bool(self.pipelines):
            self.populate_pipelines()

    def _get_email_from_config(self):
        config = self.call_carrot("carrot_cli config get")
        return config["email"]

    def _get_online_path(self, resource):
        """
        Currently carrot requires resources such as
        test and evaluation WDLs to be publicly
        available on a public repository, or at
        least be accessible to carrot's service account.
        """
        base_uri = self.repository.replace(
            "github.com", "raw.githubusercontent.com")
        uri = f"{base_uri}/{self.branch}/{resource}"
        try:
            urllib.request.urlopen(uri).getcode()
            return uri
        except HTTPError as e:
            raise ValueError(
                f"The resource {resource} is not accessible from {uri}. "
                f"Are all your files committed and pushed "
                f"to the given repository? {e}")

    @staticmethod
    def _get_wdl_outputs(wdl, include_supported_types_only=True):
        """
        Returns a dictionary whose keys are the output
        variables of the given WDL, and values are
        variable types.

        This method use womtool to extract outputs of
        a given WDL.
        """
        cmd = f"java -jar {WOMTOOL_PATH} outputs {wdl}"
        outputs_json = check_output(cmd, shell=True)
        outputs = json.loads(outputs_json)
        if include_supported_types_only:
            outputs = {k: v for k, v in outputs.items()
                       if v.lower() in SUPPORTED_RESULTS_TYPE}
        return outputs

    def load_pipelines(self):
        if os.path.isfile(self.pipelines_filename):
            with open(self.pipelines_filename, "r") as f:
                self.pipelines = json.load(f)
        else:
            self.persist_pipelines()

    def persist_pipelines(self):
        with open(self.pipelines_filename, "w") as f:
            json.dump(self.pipelines, f, default=lambda o: o.__dict__, indent="\t", sort_keys=True)

    def call_carrot(self, command):
        try:
            data = check_output(command, shell=True)
            return self.response_to_json(data)
        except subprocess.CalledProcessError as grepexc:
            print("error code", grepexc.returncode, grepexc.output)

    @staticmethod
    def response_to_json(response):
        # A hack to fix the
        response = response.decode("utf8").replace("\'", "\"").replace("None", "\"None\"")
        if "Encountered a connection error." in response:
            raise Exception("Please check your connection, are you connected to the VPN?\n" + response)
        if "\"status\": 500" in response:
            raise Exception(response)
        return json.loads(response)

    @staticmethod
    def get_timestamp():
        return calendar.timegm(time.gmtime())

    def populate_pipelines(self):
        pipelines = {}
        for root, dirs, _ in os.walk(self.working_dir, topdown=True):
            if root == self.working_dir:
                continue
            pipelines[os.path.basename(root)] = dirs.copy()
            dirs.clear()

        for pipeline, templates in pipelines.items():
            test_wdl = self._get_online_path(f"{WDLS_DIR}/{pipeline}.wdl")
            test_workflow_outputs = self._get_wdl_outputs(f"{WDLS_DIR_RELATIVE}/{pipeline}.wdl")
            if not bool(test_workflow_outputs):
                # TODO: skip the pipeline instead
                raise Exception("No supported result type for the pipeline.")

            test_results = []
            for var_name, var_type in test_workflow_outputs.items():
                result = self.create_result(name=var_name, result_type=var_type)
                result.var_name = var_name
                test_results.append(result)

            created_pipeline = self.create_pipeline()
            for template in templates:
                eval_wdl = self._get_online_path(f"{WDLS_TEST_DIR}/{pipeline}/{template}/eval.wdl")
                eval_workflow_outputs = self._get_wdl_outputs(f"{pipeline}/{template}/eval.wdl")
                if not bool(eval_workflow_outputs):
                    # TODO: skip the pipeline instead
                    raise Exception("No supported result type for the pipeline.")

                eval_results = []
                for var_name, var_type in eval_workflow_outputs.items():
                    result = self.create_result(name=var_name, result_type=var_type)
                    result.var_name = var_name
                    eval_results.append(result)

                created_template = self.create_template(created_pipeline.uuid, test_wdl, eval_wdl)
                for result in test_results + eval_results:
                    self.create_template_to_result_mapping(created_template, result)
                    created_template.results.append(result)

                created_template.test = self.create_test(created_template, os.path.join(self.working_dir, pipeline, template))
                created_pipeline.templates[created_template.uuid] = created_template

            self.pipelines[created_pipeline.uuid] = created_pipeline
        self.persist_pipelines()

    def create_pipeline(self, name=None, description=None, timestamp=True):
        name = name or "gatk_sv"
        # Name is a unique identifier in carrot,
        # if not timestamped, make sure the given
        # name does not already exist.
        if timestamp:
            name = f"{name}_{get_timestamp()}"
        description = description or "A pipeline created by the GATK-SV carrot_helper.py"
        cmd = f"carrot_cli pipeline create " \
              f"--name '{name}' " \
              f"--description '{description}'"
        response = self.call_carrot(cmd)
        return Pipeline(**response)

    def create_template(self, pipeline_id, test_wdl, eval_wdl, name=None, description=None, timestamp=True):
        name = name or "gatk_sv"
        if timestamp:
            name = f"{name}_{get_timestamp()}"
        description = description or "A template created by the GATK-SV carrot_helper.py"
        cmd = f"carrot_cli template create " \
              f"--pipeline_id {pipeline_id} " \
              f"--name {name} " \
              f"--description '{description}' " \
              f"--test_wdl {test_wdl} " \
              f"--eval_wdl {eval_wdl}"
        response = self.call_carrot(cmd)
        return Template(**response)

    def create_result(self, result_type, name=None, description=None, timestamp=True):
        """
        Outputs of test and evaluation WDLs should be
        defined as `result`s in carrot.
        """
        name = name or "gatk"
        if timestamp:
            name = f"{name}_{get_timestamp()}"
        description = description or "A results type created by GATK-sv carrot_helper.py"
        cmd = f"carrot_cli result create " \
              f"--name {name} " \
              f"--description '{description}' " \
              f"--result_type {result_type.lower()}"  # lower type is required by carrot.
        response = self.call_carrot(cmd)
        return Result(**response)

    def create_template_to_result_mapping(self, template, result):
        cmd = f"carrot_cli template map_to_result " \
              f"{template.uuid} " \
              f"{result.uuid} " \
              f"{result.var_name}"
        self.call_carrot(cmd)

    def create_test(self, template, template_path, name=None, description=None, timestamp=True):
        name = name or "gatk"
        if timestamp:
            name = f"{name}_{get_timestamp()}"
        description = description or "A test type created by GATK-sv carrot_helper.py"

        def_eval = os.path.join(template_path, DEFAULT_EVAL_INPUTS_FILENAME)
        if not os.path.isfile(def_eval):
            raise FileNotFoundError(f"Default inputs for evaluation WDL does not exist; expected file: {def_eval}")
        def_test = os.path.join(template_path, DEFAULT_TEST_INPUTS_FILENAME)
        if not os.path.isfile(def_test):
            raise FileNotFoundError(f"Default inputs for test WDL does not exist; expected file: {def_test}")

        cmd = f"carrot_cli test create" \
              f" --name '{name}' " \
              f"--description '{description}' " \
              f"--template_id {template.uuid} " \
              f"--eval_input_defaults {def_eval} " \
              f"--test_input_defaults {def_test}"
        response = self.call_carrot(cmd)
        return Test(**response, template_path=template_path)


class BaseModel:
    def __init__(self, pipeline_id, name, description, created_at, created_by):
        self.uuid = pipeline_id
        self.name = name
        self.description = description
        self.created_at = created_at
        self.created_by = created_by


class Pipeline(BaseModel):
    """
    Carrot pipeline groups tests together.
    """
    def __init__(self, pipeline_id, name, description, created_at, created_by):
        super().__init__(pipeline_id, name, description, created_at, created_by)
        self.templates = {}


class Template(BaseModel):
    """
    A carrot template is composed of WDL to be tested,
    a WDL that evaluates the tested WDL's output, and
    a set of inputs for both test and evaluation WDLs.
    """
    def __init__(self, template_id, pipeline_id, test_wdl, eval_wdl, name, description, created_at, created_by):
        super().__init__(template_id, name, description, created_at, created_by)
        self.pipeline_id = pipeline_id
        self.test_wdl = test_wdl
        self.eval_wdl = eval_wdl
        self.results = []
        self.test = None


class Result(BaseModel):
    def __init__(self, result_id, result_type, name, description, created_at, created_by):
        super().__init__(result_id, name, description, created_at, created_by)
        self.var_name = None
        self.result_type = result_type


class Test(BaseModel):
    def __init__(self, test_id, name, template_id, template_path, eval_input_defaults, test_input_defaults, description, created_at, created_by):
        super().__init__(test_id, name, description, created_at, created_by)
        self.template_id = template_id
        self.eval_input_defaults = eval_input_defaults
        self.test_input_defaults = test_input_defaults
        self.template_path = template_path


def get_timestamp():
    return calendar.timegm(time.gmtime())


def response_to_json(response):
    return json.loads(response.decode("utf8").replace("\'", "\"").replace("None", "\"None\""))


def call_carrot_cli(command):
    # process = Popen(command.split(" "), stdout=subprocess.PIPE)
    # data, err = process.communicate()
    try:
        data = check_output(command, shell=True)
        # return response_to_json(data) if process.returncode == 0 else err
        return response_to_json(data)
    except subprocess.CalledProcessError as grepexc:
        print("error code", grepexc.returncode, grepexc.output)


def create_test():
    pipeline_id = "2214f31b-8aee-48cd-a81f-a0bd3091fa54"
    test_wdl = "https://raw.githubusercontent.com/VJalili/gatk-sv/ehdn_unit_test/wdl/ExpansionHunterDenovo.wdl"
    eval_wdl = "https://raw.githubusercontent.com/VJalili/gatk-sv/ehdn_unit_test/wdl_test/ExpansionHunterDenovo/eval_casecontrol_locus.wdl"
    report_id = "ae61fbff-4270-4528-91a4-1c0091076289"
    test_input_defaults = "wdl_test/ExpansionHunterDenovo/test_casecontrol_locus_default_inputs.json"
    eval_input_defaults = "wdl_test/ExpansionHunterDenovo/eval_casecontrol_locus_default_inputs.json"
    results = {
        "EHdnSTRAnalysis.multisample_profile": "b4a2bdf7-a6b1-479d-82a5-53f63ac9861c",
        "EvalCaseControlLocus.data_file": "3ecb2895-50a0-40d0-8931-6d8d6cfd2b5c"
    }

    cmd = f"carrot_cli template create --pipeline_id {pipeline_id} --name test_{get_timestamp()} --test_wdl {test_wdl} --eval_wdl {eval_wdl}"
    output = call_carrot_cli(cmd)
    template_id = output["template_id"]
    print(f"\nCreated Template:\n{output}\n")

    for workflow_output_label, result_id in results.items():
        cmd = f"carrot_cli template map_to_result {template_id} {result_id} {workflow_output_label}"
        output = call_carrot_cli(cmd)
        print(f"\nCreated Result:\n{output}\n")

    cmd = f"carrot_cli template map_to_report {template_id} {report_id}"
    output = call_carrot_cli(cmd)
    print(f"\nCreated Template mapping to Report:\n{output}\n")

    cmd = f"carrot_cli test create --name test_{get_timestamp()} --template_id {template_id} --test_input_defaults {test_input_defaults} --eval_input_defaults {eval_input_defaults}"
    output = call_carrot_cli(cmd)
    print(f"\nCreated Test:\n{output}\n")
    return output["test_id"]


def run_test(test_id):
    test_input = "wdl_test/ExpansionHunterDenovo/test_casecontrol_locus_inputs.json"
    eval_input = "wdl_test/ExpansionHunterDenovo/eval_casecontrol_locus_inputs.json"
    cmd = f"carrot_cli test run --test_input {test_input} --eval_input {eval_input} {test_id}"
    output = call_carrot_cli(cmd)
    print(f"\nExecuted Test:\n{output}\n")
    print(f"\n\nRun ID: {output['run_id']}")
    print(f"Get Run Status: carrot_cli run find_by_id {output['run_id']}")


def delete_all_software_created_by(email="jalili.vahid@broadinstitute.org"):
    command = f"carrot_cli software find --created_by {email}"
    response = check_output(command, shell=True)
    my_softwares = json.loads(response.decode("utf8").replace("\'", "\"").replace("None", "\"None\""))
    for x in my_softwares:
        response = check_output(command, shell=True)


if __name__ == '__main__':

    carrot_helper = CarrotHelper()
    exit()

    parent_parser = argparse.ArgumentParser(description="Helper methods to scaffold and updated carrot tests.")
    subparsers = parent_parser.add_subparsers(title="service", dest="service_command")

    update_parser = subparsers.add_parser("update-and-run", help="Update test")
    # update_parser.add_argument("target_wdl", help="The WDL whose test should be updated.")

    args = parent_parser.parse_args()
    test_id = create_test()
    run_test(test_id)
