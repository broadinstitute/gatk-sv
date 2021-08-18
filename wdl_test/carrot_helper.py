"""
Requirements:
- Java
- Womtool
- carrot_cli
"""

import argparse
import hashlib
import subprocess
import sys
from subprocess import check_output
import time
import calendar
import json
import os
import shutil
import urllib.request
from dataclasses import dataclass, field
from typing import List
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

# A JSON file to which the current setup is
# persisted-to and loaded-from.
PIPELINES_JSON = ".carrot_pipelines.json"

# This file contains all the submitted runs.
# If renamed, also rename in `.gitignore`.
RUNS_JSON = ".runs.json"

# This is the list of the currently
# supported results type by Carrot.
# Keep the list lowercase, as carrot is
# case-sensitive and requires lower-case. See:
# https://github.com/broadinstitute/carrot/blob/1b2905e634e0887a94c63743a78b1731bd3a637c/src/custom_sql_types.rs#L158-L162
SUPPORTED_RESULTS_TYPE = ["numeric", "file", "text"]


class CarrotHelper:
    def __init__(self, working_dir=".",
                 pipelines_filename=PIPELINES_JSON,
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
        else:
            self.update_resources()

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
    def _get_test_wdl_local_path(wdl):
        return f"{WDLS_DIR_RELATIVE}/{wdl}.wdl"

    @staticmethod
    def _get_eval_wdl_local_path(pipeline, template):
        return f"{pipeline}/{template}/eval.wdl"

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

    @staticmethod
    def get_checksum(filename):
        """
        Implemented based on: https://stackoverflow.com/a/22058673/947889
        """
        buffer_size = 65536  # 64kb
        sha3 = hashlib.sha3_256()
        with open(filename, "rb") as f:
            while True:
                data = f.read(buffer_size)
                if not data:
                    break
                sha3.update(data)
        return sha3.hexdigest()

    def load_pipelines(self):
        if os.path.isfile(self.pipelines_filename):
            with open(self.pipelines_filename, "r") as f:
                self.pipelines = json.loads(f.read())
                # At this point, pipelines is loaded as a
                # dictionary of dictionaries, and the
                # following code converts it the object
                # types as it was when serialized.
                for k, v in self.pipelines.items():
                    pipeline = Pipeline(**v)
                    for tk, tv in pipeline.templates.items():
                        template = Template(**tv)
                        results = []
                        for result in template.results:
                            results.append(Result(**result))
                        template.results = results
                        template.test = Test(**template.test)
                        pipeline.templates[tk] = template
                    self.pipelines[k] = pipeline
        else:
            self.persist_object(self.pipelines, self.pipelines_filename)

    def update_resources(self):
        for pipeline_key, pipeline_value in self.pipelines.items():
            self.setup_pipeline(pipeline_key, pipeline_value, pipeline_value.templates.keys())
        self.persist_object(self.pipelines, self.pipelines_filename)

    def setup_pipeline(self, pipeline_key, pipeline_value, templates):
        test_wdl_local = self._get_test_wdl_local_path(pipeline_key)
        test_workflow_outputs = self._get_wdl_outputs(test_wdl_local)
        if not bool(test_workflow_outputs):
            # TODO: skip the pipeline instead
            raise Exception("No supported result type for the pipeline.")

        if not pipeline_value:
            pipeline_value = self.create_pipeline()
        test_wdl = self._get_online_path(f"{WDLS_DIR}/{pipeline_key}.wdl")
        test_wdl_checksum = self.get_checksum(test_wdl_local)
        for template in templates:
            if template in pipeline_value.templates and (pipeline_value.templates[template].test_wdl_checksum == test_wdl_checksum and pipeline_value.templates[template].eval_wdl_checksum == self.get_checksum(self._get_eval_wdl_local_path(pipeline_key, template))):
                continue
            created_template = self.setup_template(pipeline_key, pipeline_value, template, test_wdl, test_wdl_checksum, test_workflow_outputs) #, test_results)
            pipeline_value.templates[template] = created_template

        return pipeline_value

    def get_pipelines(self, dirs):
        """
        Takes a list of directories as input, traverses every
        directory and sub-directories and extracts pipelines.
        """
        pipelines = {}
        wd = self.working_dir
        dirs = [x for x in dirs if os.path.isdir(x)]
        for d in dirs:
            d = os.path.normpath(d).split(os.sep)
            if len(d) == 0:
                continue
            pipeline = d[0]
            templates = [d[1]] if len(d) == 2 else None
            run_dirs = [d[2]] if len(d) == 3 else None

            if pipeline not in pipelines:
                pipelines[pipeline] = {}
            templates = templates or [x for x in next(os.walk(os.path.join(wd, pipeline)))[1]]
            for t in templates:
                if t not in pipelines[pipeline]:
                    pipelines[pipeline][t] = []
                run_dirs = run_dirs or [x for x in next(os.walk(os.path.join(wd, pipeline, t)))[1]]
                for r in run_dirs:
                    if r not in pipelines[pipeline][t]:
                        pipelines[pipeline][t].append(r)
        return pipelines

    @staticmethod
    def persist_object(obj, filename):
        with open(filename, "w") as f:
            json.dump(obj, f, default=lambda o: o.__dict__, indent="\t", sort_keys=True)

    def call_carrot(self, command, resource_id=None):
        """

        :param command:
        :param resource_id: If provided, the provided ID field will be renamed to `uuid`.
        For instance; `resource_id="pipeline_id"` will lead to changing `pipeline_id` to `uuid`
        in the deserialized object.
        :return:
        """
        try:
            data = check_output(command, shell=True)
            return self.response_to_json(data, resource_id)
        except subprocess.CalledProcessError as grepexc:
            print("error code", grepexc.returncode, grepexc.output)

    @staticmethod
    def response_to_json(response, resource_id):
        # A hack to fix the
        response = response.decode("utf8").replace("\'", "\"").replace("None", "\"None\"")
        if "Encountered a connection error." in response:
            raise Exception("Please check your connection, are you connected to the VPN?\n" + response)
        if "\"status\": 500" in response:
            raise Exception(response)
        j = json.loads(response)
        if resource_id:
            j["uuid"] = j.pop(resource_id)
        return j

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
            self.pipelines[pipeline] = self.setup_pipeline(pipeline, None, templates)
        self.persist_object(self.pipelines, self.pipelines_filename)

    def create_pipeline(self, name=None, description=None, timestamp=True):
        name = name or "gatk_sv"
        # Name is a unique identifier in carrot,
        # if not timestamped, make sure the given
        # name does not already exist.
        if timestamp:
            name = f"{name}_{get_timestamp()}"
        description = description or "A pipeline created by the GATK-SV carrot_helper.py"
        cmd = f"carrot_cli -q pipeline create " \
              f"--name '{name}' " \
              f"--description '{description}'"
        response = self.call_carrot(cmd, "pipeline_id")
        return Pipeline(**response)

    def setup_results(self, results, workflow_outputs):
        updated_results = []
        for var_name, var_type in workflow_outputs.items():
            # because carrot lowers output types.
            var_type = var_type.lower()
            for result in results:
                if result.var_name == var_name and result.result_type == var_type:
                    updated_results.append(result)
                    break
            else:
                result = self.create_result(name=var_name, result_type=var_type)
                result.var_name = var_name
                updated_results.append(result)
        return updated_results

    def setup_template(self, pipeline, created_pipeline, template, test_wdl, test_wdl_checksum, test_workflow_outputs): #, test_results):
        eval_wdl = self._get_online_path(f"{WDLS_TEST_DIR}/{pipeline}/{template}/eval.wdl")
        eval_wdl_local = self._get_eval_wdl_local_path(pipeline, template)
        eval_wdl_checksum = self.get_checksum(eval_wdl_local)
        eval_workflow_outputs = self._get_wdl_outputs(eval_wdl_local)
        if not bool(eval_workflow_outputs):
            # TODO: skip the pipeline instead
            raise Exception("No supported result type for the pipeline.")

        created_template = self.create_template(created_pipeline.uuid, test_wdl, test_wdl_checksum, eval_wdl, eval_wdl_checksum)
        test_results = self.setup_results(created_pipeline.templates[template].results, test_workflow_outputs)
        eval_results = self.setup_results(created_pipeline.templates[template].results, eval_workflow_outputs)
        for result in test_results + eval_results:
            self.create_template_to_result_mapping(created_template, result)
            created_template.results.append(result)

        created_template.test = self.create_test(created_template, os.path.join(self.working_dir, pipeline, template))
        return created_template

    def create_template(self, pipeline_id, test_wdl, test_wdl_checksum, eval_wdl, eval_wdl_checksum, name=None, description=None, timestamp=True):
        name = name or "gatk_sv"
        if timestamp:
            name = f"{name}_{get_timestamp()}"
        description = description or "A template created by the GATK-SV carrot_helper.py"
        cmd = f"carrot_cli -q template create " \
              f"--pipeline_id {pipeline_id} " \
              f"--name {name} " \
              f"--description '{description}' " \
              f"--test_wdl {test_wdl} " \
              f"--eval_wdl {eval_wdl}"
        response = self.call_carrot(cmd, "template_id")
        return Template(**response, test_wdl_checksum=test_wdl_checksum, eval_wdl_checksum=eval_wdl_checksum)

    def create_result(self, result_type, name=None, description=None, timestamp=True):
        """
        Outputs of test and evaluation WDLs should be
        defined as `result`s in carrot.
        """
        name = name or "gatk"
        if timestamp:
            name = f"{name}_{get_timestamp()}"
        description = description or "A results type created by GATK-sv carrot_helper.py"
        cmd = f"carrot_cli -q result create " \
              f"--name {name} " \
              f"--description '{description}' " \
              f"--result_type {result_type.lower()}"  # lower type is required by carrot.
        response = self.call_carrot(cmd, "result_id")
        return Result(**response)

    def create_template_to_result_mapping(self, template, result):
        cmd = f"carrot_cli -q template map_to_result " \
              f"{template.uuid} " \
              f"{result.uuid} " \
              f"{result.var_name}"
        self.call_carrot(cmd)

    def create_test(self, template, template_path, name=None, description=None, timestamp=True):
        def_eval = os.path.join(template_path, DEFAULT_EVAL_INPUTS_FILENAME)
        if not os.path.isfile(def_eval):
            raise FileNotFoundError(f"Default inputs for evaluation WDL does not exist; expected file: {def_eval}")
        def_test = os.path.join(template_path, DEFAULT_TEST_INPUTS_FILENAME)
        if not os.path.isfile(def_test):
            raise FileNotFoundError(f"Default inputs for test WDL does not exist; expected file: {def_test}")
        def_eval_checksum = self.get_checksum(def_eval)
        def_test_checksum = self.get_checksum(def_test)

        if template.test.test_input_defaults_checksum == def_test_checksum and template.test.eval_input_defaults_checksum == def_eval_checksum:
            return template.test

        name = name or "gatk"
        if timestamp:
            name = f"{name}_{get_timestamp()}"
        description = description or "A test type created by GATK-sv carrot_helper.py"
        cmd = f"carrot_cli -q test create " \
              f"--name '{name}' " \
              f"--description '{description}' " \
              f"--template_id {template.uuid} " \
              f"--eval_input_defaults {def_eval} " \
              f"--test_input_defaults {def_test}"
        response = self.call_carrot(cmd, "test_id")
        return Test(**response, template_path=template_path, eval_input_defaults_checksum=def_eval_checksum, test_input_defaults_checksum=def_test_checksum)

    def create_run(self, test, run_dir, name=None, timestamp=True):
        name = name or "gatk"
        if timestamp:
            name = f"{name}_{get_timestamp()}"
        eval_input = os.path.join(test.template_path, run_dir, EVAL_INPUTS_FILENAME)
        test_input = os.path.join(test.template_path, run_dir, TEST_INPUTS_FILENAME)
        cmd = f"carrot_cli -q test run --name {name} --test_input {test_input} --eval_input {eval_input} {test.uuid}"
        response = self.call_carrot(cmd)
        return Run(**response, run_dir=os.path.join(test.template_path, run_dir))

    def update_runs_status(self):
        try:
            with open(RUNS_JSON, "r") as runs_json:
                runs_dict = json.load(runs_json)
                runs = {}
                for k, v in runs_dict.items():
                    runs[k] = Run(**v)
        # TODO: instead of capturing exceptions and printing messages here
        # a better alternative would be to capture the errors at the callers level.
        except FileNotFoundError:
            print(f"The runs JSON file, {RUNS_JSON}, does not exist.")
        except json.decoder.JSONDecodeError:
            print(f"The runs JSON file, {RUNS_JSON}, is empty or "
                  f"corrupted; have you submitted any test runs?")

        updated_status = []
        for run_id, run in runs.items():
            if run.status in ["succeeded", "failed"]:
                updated_status.append(RunStatus(run_id, run.run_dir, run.status, run.status, run.test_cromwell_job_id, run.eval_cromwell_job_id))
                continue
            cmd = f"carrot_cli run find_by_id {run_id}"
            previous_status = run.status
            # At time of submitting a run, the run may not
            # have a test cromwell job ID and certainly not
            # eval cromwell job ID, but updated status, depending
            # on the status of test execution, can update both
            # of these IDs, which are useful for debugging
            # failed tests.
            updated_run = self.call_carrot(cmd)
            run.__dict__.update(updated_run)
            updated_status.append(RunStatus(run_id, run.run_dir, previous_status, run.status, run.test_cromwell_job_id or "None", run.eval_cromwell_job_id or "None"))
            runs[run_id] = run
        self.persist_object(runs, RUNS_JSON)
        return updated_status

    def submit_tests(self, tests_dir):
        try:
            runs = {}
            with open(RUNS_JSON, "r") as runs_json:
                runs = json.load(runs_json)
        except json.decoder.JSONDecodeError:
            # Happens when this script is run for the
            # first time and the file does not exist.
            # Passing this exception should be fine
            # since a JSON object will be created and
            # persisted to the file if creating a `run`
            # was run successfully.
            pass

        c = 0
        print("Submitting tests ...")
        print('%-40s%-80s' % ("Run ID", "Test"))
        pipelines = self.get_pipelines(tests_dir)
        for pipeline, templates in pipelines.items():
            for template, run_dirs in templates.items():
                for run_dir in run_dirs:
                    c += 1
                    test = self.pipelines[pipeline].templates[template].test
                    run = self.create_run(test, run_dir)
                    runs[run.uuid] = run
                    self.persist_object(runs, RUNS_JSON)
                    print('%-40s%-80s' % (run.uuid, pipeline + "/" + template + "/" + run_dir))
        print(f"Submitted {c} test(s).")

    @staticmethod
    def pretty_print_runs_status(statuses):
        def format_run_dir(txt):
            if run_dir_col_max_width < 3:
                return ""
            if run_dir_col_max_width < 10:
                return "..."
            return "..." + txt[-(run_dir_col_max_width - 3):] if len(txt) > run_dir_col_max_width else txt

        t_size = shutil.get_terminal_size()
        # 36: UUID columns
        # 16: status columns
        # 15: (5 columns - 1) * 3, spaces between columns
        run_dir_col_max_width = t_size.columns - (((36 * 3) + (16 * 2)) + 15)
        row_format = f"{{:<36}}   {{:<{run_dir_col_max_width}}}   {{:<16}}   {{:<16}}   {{:<36}}   {{:<36}}"

        print(row_format.format(
            "Run UUID",
            "Run Directory",
            "Previous Status",
            "Latest Status",
            "Test Cromwell Job ID",
            "Eval Cromwell Job ID"))

        for status in statuses:
            ec = "\033[0m"
            if status.latest_status == "succeeded":
                sc = "\033[92m"  # green color
            elif status.latest_status == "failed":
                sc = "\033[91m"  # Red color
            else:
                sc = "\033[93m"  # yellow color

            print(row_format.format(
                status.run_uuid,
                format_run_dir(status.run_dir),
                status.previous_status,
                f"{sc}{status.latest_status}{ec}",
                status.test_cromwell_job_id,
                status.eval_cromwell_job_id))


@dataclass
class BaseModel:
    uuid: str
    name: str
    created_at: str
    created_by: str
    description: str


@dataclass
class Pipeline(BaseModel):
    """
    Carrot pipeline groups tests together.
    """
    # See the following docs on mutable default values:
    # https://docs.python.org/3/library/dataclasses.html#mutable-default-values
    templates: dict = field(default_factory=dict)


@dataclass
class Template(BaseModel):
    """
    A carrot template is composed of WDL to be tested,
    a WDL that evaluates the tested WDL's output, and
    a set of inputs for both test and evaluation WDLs.
    """
    pipeline_id: str
    test_wdl: str
    test_wdl_checksum: str
    eval_wdl: str
    eval_wdl_checksum: str
    results: List = field(default_factory=list)
    test: str = None


@dataclass
class Result(BaseModel):
    result_type: str
    var_name: str = None


@dataclass
class Test(BaseModel):
    template_id: str
    eval_input_defaults_checksum: str
    test_input_defaults_checksum: str
    template_path: str
    eval_input_defaults: dict = field(default_factory=dict)
    test_input_defaults: dict = field(default_factory=dict)


class Run(BaseModel):
    def __init__(self, name, test_id, test_cromwell_job_id,
                 eval_cromwell_job_id, status, finished_at,
                 test_input, eval_input, created_at, created_by,
                 description=None, run_id=None, uuid=None,
                 results=None, run_dir=None):
        super().__init__(run_id or uuid, name, description,
                         created_at, created_by)
        self.test_id = test_id
        self.test_cromwell_job_id = test_cromwell_job_id
        self.eval_cromwell_job_id = eval_cromwell_job_id
        self.status = status
        self.finished_at = finished_at
        self.test_input = test_input
        self.eval_input = eval_input
        self.results = results
        self.run_dir = run_dir


@dataclass
class RunStatus:
    run_uuid: str
    run_dir: str
    previous_status: str
    latest_status: str
    test_cromwell_job_id: str
    eval_cromwell_job_id: str


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


def main():
    wd = "."
    parser = argparse.ArgumentParser(
        prog=sys.argv[0],
        description="Helper methods to create, "
                    "setup, and run carrot tests.")
    subparsers = parser.add_subparsers(help='Sub-commands')

    test_parser = subparsers.add_parser("test", help="Create, list, and run tests.")
    test_subparsers = test_parser.add_subparsers(title="sub-commands for test", dest="test", help="Commands to interact with carrot tests.")
    test_run_parser = test_subparsers.add_parser("run")
    test_run_parser.add_argument("tests_dir", nargs="+")
    test_subparsers.add_parser("update_status")

    args = parser.parse_args()

    carrot_helper = CarrotHelper(working_dir=wd)
    if "test" in args:
        args_dict = vars(args)
        if args.test == "run":
            carrot_helper.submit_tests(args_dict.get("tests_dir"))
        elif args.test == "update_status":
            updated_status = carrot_helper.update_runs_status()
            carrot_helper.pretty_print_runs_status(updated_status)


if __name__ == '__main__':
    main()
