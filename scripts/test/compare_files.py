import argparse
import gzip
import json
import os
from itertools import chain
from subprocess import DEVNULL, STDOUT, check_call


class BaseCompareAgent:
    def __init__(self, working_dir=None):
        self.working_dir = working_dir

    def get_filename(self, obj):
        return obj.replace("gs://", self.working_dir)

    def get_obj(self, obj):
        """
        Ensures the given Google Cloud Storage
        object (obj) is available in the working
        directory, and returns its filename in
        the working directory.
        """
        raise NotImplementedError


class VCFCompareAgent(BaseCompareAgent):
    def __init__(self, working_dir=None):
        super(VCFCompareAgent, self).__init__(working_dir)

        # Delimiter
        self.d = "\t"
        self.id_col = 2

    def get_obj(self, obj):
        """
        Ensures the given VCF object is
        available in the working directory:
        if it exists, returns its filename, and
        If it does not, downloads the VCF object and
        returns its filename.
        """
        filename = self.get_filename(obj)
        if not os.path.isfile(filename):
            if not os.path.isfile(filename):
                check_call(
                    ["gsutil", "-m", "cp", obj, filename],
                    stdout=DEVNULL, stderr=STDOUT)
        return filename

    def equals(self, x, y):
        """
        Gets two VCF objects (Google Cloud Storage URI),
        x and y, and returns true if files are identical,
        and false if otherwise. Additionally, it returns the
        compared files.
        """
        x = self.get_obj(x)
        y = self.get_obj(y)

        with gzip.open(x, "rt", encoding="utf-8") as X, \
             gzip.open(y, "rt", encoding="utf-8") as Y:
            # Not an optimal search, but it works
            # if the files are not sorted, which
            # is expected for some files in the
            # pipeline.
            for x_line in X:
                if x_line.startswith("chr"):
                    x_columns = x_line.split(self.d)
                    for y_line in Y:
                        if y_line.startswith("chr"):
                            y_columns = y_line.split(self.d)
                            if len(x_columns) == len(y_columns):
                                for c in chain(range(0, self.id_col), range(self.id_col + 1, len(x_columns))):
                                    if x_columns[c] != y_columns[c]:
                                        break
                                else:
                                    break
                    else:
                        return False, x, y
        return True, x, y


class CompareWorkflowOutputs:
    def __init__(self, working_dir):
        self.working_dir = working_dir
        self.filetypes_to_compare = {
            "vcf.gz": VCFCompareAgent(self.working_dir)
        }

    def _get_filtered_outputs(self, task_outputs):
        """
        Iterates through the outputs of a task and
        filters the outputs whose file type match
        the types subject to comparison (i.e.,
        types defined in filetypes_to_compare).
        """
        filtered_outputs = {}
        if not isinstance(task_outputs, list):
            task_outputs = [task_outputs]

        for task_output in task_outputs:
            if not isinstance(task_output, str):
                # Happens when output is not a file,
                # e.g., when it is a number.
                continue
            for ext in self.filetypes_to_compare:
                if task_output.endswith(ext):
                    filtered_outputs[ext] = task_output
        return filtered_outputs

    def _traverse_workflow(
            self, calls, parent_workflow="",
            traverse_sub_workflows=False):
        output_files = {}

        def update_output_files(outputs):
            if len(outputs) > 0:
                output_files[
                    ((parent_workflow + ".") if parent_workflow else "") +
                    f"{workflow}.{out_label}" +
                    (("." + str(run["shardIndex"])) if run["shardIndex"] != -1 else "")] \
                    = outputs

        for workflow, runs in calls.items():
            for run in runs:
                if run["executionStatus"] != "Done":
                    continue
                for out_label, out_files in run["outputs"].items():
                    if not out_files:
                        continue
                    update_output_files(self._get_filtered_outputs(out_files))
                if traverse_sub_workflows and "subWorkflowMetadata" in run:
                    update_output_files(self._get_filtered_outputs(out_files))
                    output_files.update(
                        self._traverse_workflow(
                            run["subWorkflowMetadata"]["calls"],
                            workflow, traverse_sub_workflows))
        return output_files

    def _get_output_files(self, filename, traverse_sub_workflows=False):
        """
        Iterates through a given cromwell metadata file
        and filters the output files to be compared.
        """
        with open(filename, "r") as f:
            metadata = json.load(f)
            output_files = self._traverse_workflow(
                metadata["calls"],
                traverse_sub_workflows=traverse_sub_workflows)
        return output_files

    def get_mismatches(
            self, reference_metadata, target_metadata,
            traverse_sub_workflows=False):
        """
        Takes two metadata files (both belonging to a common
        workflow execution), iterates through the outputs of
        their task, downloads the objects if not already exist
        in the working directory, compares the corresponding
        files, and returns the files that do not match.
        """
        ref_output_files = self._get_output_files(reference_metadata, traverse_sub_workflows)
        test_output_files = self._get_output_files(target_metadata, traverse_sub_workflows)

        # For coloring the prints; see the following SO
        # answer for details: https://stackoverflow.com/a/287944/947889
        color_green = '\033[92m'
        color_red = '\033[91m'
        color_endc = '\033[0m'

        mismatches = {}
        i = 0
        for call, ref_outputs in ref_output_files.items():
            i += 1
            print(f"Comparing\t{i}/{len(ref_output_files)}\t{call} ... ", end="")
            if len(ref_outputs) > 1:
                raise NotImplementedError
            for extension, obj in ref_outputs.items():
                equals, x, y = \
                    self.filetypes_to_compare[extension].equals(
                        obj, test_output_files[call][extension])
                if not equals:
                    if call not in mismatches:
                        mismatches[call] = []
                    mismatches[call].append([x, y])
                    print(f"{color_red}mismatch{color_endc}")
                else:
                    print(f"{color_green}match{color_endc}")
        return mismatches


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Takes two cromwell metadata files as input, "
                    "reference and target, compares their corresponding "
                    "output files, and reports the files that do not match. "
                    "The two metadata files should belong to the execution "
                    "of a common workflow (e.g., one workflow with different "
                    "inputs). The script requires `gsutil` and `gzip` to be "
                    "installed and configured."
                    "\n\n"
                    "The currently supported file types are as follows."
                    "\n\t- VCF (.vcf.gz): The non-header lines of VCF files"
                    "are compared; except for the ID column, all the other "
                    "columns of a variation are expected to be identical. "
                    "The two files are expected to be equally ordered (i.e., "
                    "n-th variation in one file is compared to the "
                    "n-th variation on the other file).",
        formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument(
        "reference_metadata",
        help="Reference cromwell metadata file.")
    parser.add_argument(
        "target_metadata",
        help="Target cromwell metadata file.")
    parser.add_argument(
        "-w", "--working_dir",
        help="The directory where the files will "
             "be downloaded; default is the "
             "invocation directory.")
    parser.add_argument(
        "-o", "--output",
        help="Output file to store mismatches "
             "(in JSON format); defaults to `output.json`.")
    parser.add_argument(
        "-d", "--deep",
        action="store_true",
        help="Include sub-workflows traversing the metadata files.")

    args = parser.parse_args()

    wd = args.working_dir if args.working_dir else "."
    comparer = CompareWorkflowOutputs(wd)
    mismatches = comparer.get_mismatches(
        args.reference_metadata,
        args.target_metadata,
        args.deep)

    print(f"{len(mismatches)} files did not match.")

    output_file = \
        args.output if args.output else \
        os.path.join(wd, "output.json")
    with open(output_file, "w") as f:
        json.dump(mismatches, f, indent=2)

    print(f"Mismatches are persisted in {output_file}.")
