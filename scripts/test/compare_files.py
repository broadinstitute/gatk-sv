import argparse
import gzip
import json
import os
from metadata import ITaskOutputFilters, Metadata
from subprocess import DEVNULL, STDOUT, check_call


# For coloring the prints; see the following SO
# answer for details: https://stackoverflow.com/a/287944/947889
COLOR_ENDC = "\033[0m"
COLOR_ULINE = "\033[04m"
COLOR_BLINKING = "\033[05m"
COLOR_RED = "\033[91m"
COLOR_GREEN = "\033[92m"
COLOR_YELLOW = "\033[93m"


class FilterBasedOnExtensions(ITaskOutputFilters):

    def __init__(self, extensions):
        self.extensions = extensions

    def filter(self, metadata, outputs):
        """
        Iterates through the outputs of a task and
        filters the outputs whose file type match
        the types subject to comparison (i.e.,
        types defined in filetypes_to_compare).

        :return: An array of the filtered outputs.
        """
        filtered_outputs = {}
        if not isinstance(outputs, list):
            outputs = [outputs]

        for task_output in outputs:
            if not isinstance(task_output, str):
                # Happens when output is not a file,
                # e.g., when it is a number.
                continue
            for ext in self.extensions:
                if task_output.endswith(ext):
                    if ext not in filtered_outputs:
                        filtered_outputs[ext] = []
                    filtered_outputs[ext].append(task_output)
        return filtered_outputs


class BaseCompareAgent:
    def __init__(self, working_dir=None):
        self.working_dir = working_dir

    def get_filename(self, obj):
        return obj.replace("gs://", os.path.join(self.working_dir, ""))

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
            for x_line, y_line in zip(X, Y):
                if x_line.startswith("#") and y_line.startswith("#"):
                    continue

                x_columns = x_line.strip().split(self.d)
                y_columns = y_line.strip().split(self.d)

                if len(x_columns) != len(y_columns):
                    return False, x, y

                if any(x_columns[c] != y_columns[c]
                       for c in range(0, len(x_columns))
                       if c != self.id_col):
                    return False, x, y
        return True, x, y


class CompareWorkflowOutputs:
    def __init__(self, working_dir):
        self.working_dir = working_dir
        self.filetypes_to_compare = {
            "vcf.gz": VCFCompareAgent(self.working_dir)
        }

    def get_mismatches(self, reference_metadata,
                       target_metadata,
                       traverse_sub_workflows=False):
        """
        Takes two metadata files (both belonging to a common
        workflow execution), iterates through the outputs of
        their task, downloads the objects if not already exist
        in the working directory, compares the corresponding
        files, and returns the files that do not match.
        """
        def record_compare_result(match, reference, target):
            if not match:
                if call not in mismatches:
                    mismatches[call] = []
                mismatches[call].append([reference, target])

        # First we define a method that takes a list
        # of a task outputs, and keeps only those that
        # are files and their extension match the
        # file types that we want to compare
        # (e.g., filter only VCF files).
        filter_method = FilterBasedOnExtensions(
            self.filetypes_to_compare.keys()).filter

        # Then we create two instances of the Metadata
        # class, one for each metadata file, and we
        # invoke the `get_outputs` method which traverses
        # the outputs of task, and returns those filtered
        # by the above-defined filter.
        ref_output_files = Metadata(reference_metadata).get_outputs(
            traverse_sub_workflows, filter_method)
        test_output_files = Metadata(target_metadata).get_outputs(
            traverse_sub_workflows, filter_method)

        mismatches = {}
        i = 0

        r_t = ref_output_files.keys() - test_output_files.keys()
        t_r = test_output_files.keys() - ref_output_files.keys()
        if r_t or t_r:
            print(f"\n{COLOR_BLINKING}WARNING!{COLOR_ENDC}")
            print(f"The reference and test metadata files differ "
                  f"in their outputs; "
                  f"{COLOR_ULINE}the differences will be skipped.{COLOR_ENDC}")
            if r_t:
                print(f"\t{len(r_t)}/{len(ref_output_files.keys())} "
                      f"outputs of the reference are not in the test:")
                for x in r_t:
                    print(f"\t\t- {x}")
            if t_r:
                print(f"\t{len(t_r)}/{len(test_output_files.keys())} "
                      f"outputs of the test are not in the reference:")
                for x in t_r:
                    print(f"\t\t- {x}")
            print("\n")

        [ref_output_files.pop(x) for x in r_t]
        print(f"{COLOR_YELLOW}Comparing {len(ref_output_files)} "
              f"files that are common between reference and test "
              f"metadata files and their respective task is executed "
              f"successfully.{COLOR_ENDC}")
        for call, ref_outputs in ref_output_files.items():
            i += 1
            matched = True
            print(f"Comparing\t{i}/{len(ref_output_files)}\t{call} ... ", end="")
            for extension, objs in ref_outputs.items():
                if len(objs) != len(test_output_files[call][extension]):
                    record_compare_result(False, objs, test_output_files[call][extension])
                    matched = False
                    continue
                for idx, obj in enumerate(objs):
                    equals, x, y = \
                        self.filetypes_to_compare[extension].equals(
                            obj, test_output_files[call][extension][idx])
                    record_compare_result(equals, x, y)
                    if not equals:
                        matched = False
            if matched:
                print(f"{COLOR_GREEN}match{COLOR_ENDC}")
            else:
                print(f"{COLOR_RED}mismatch{COLOR_ENDC}")
        return mismatches


def main():
    parser = argparse.ArgumentParser(
        description="Takes two cromwell metadata files as input, "
                    "reference and target, compares their corresponding "
                    "output files, and reports the files that do not match. "
                    "The two metadata files should belong to the execution "
                    "of a common workflow (e.g., one workflow with different "
                    "inputs). The script requires `gsutil` and `gzip` to be "
                    "installed and configured. If the output of a task is an "
                    "array of files, the reference and target arrays are "
                    "expected to be in the same order."
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

    if len(mismatches) == 0:
        print(f"{COLOR_GREEN}All the compared files matched.{COLOR_ENDC}")
    else:
        print(f"{COLOR_RED}{len(mismatches)} of the compared files did not match.{COLOR_ENDC}")
        output_file = \
            args.output if args.output else \
            os.path.join(wd, "output.json")
        with open(output_file, "w") as f:
            json.dump(mismatches, f, indent=2)
        print(f"Mismatches are persisted in {output_file}.")


if __name__ == '__main__':
    main()
