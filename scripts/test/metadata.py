import json
import types
from abc import ABC, abstractmethod


class ITaskOutputFilters(ABC):
    """
    An interface that should be implemented by
    custom filtering methods to be used with Metadata.

    This design follows the principles of strategy pattern,
    where a custom method can be used to augment the default
    behavior of an algorithm. Here, this design is used to
    decouple the filtering of tasks outputs (e.g., only extract
    files with certain extension) from metadata traversal.
    """

    @abstractmethod
    def filter(self, metadata, outputs):
        """
        How to filter the output of a task.

        Note that the method is stateful; i.e.,
        it has references to both self and to
        the instance of Metadata class that
        invokes this method.

        :param metadata: `self` of the instance
        of the Metadata class that calls this method.

        :param outputs: The values of a key in the
        `outputs` field in a metadata file. e.g.,
        `metadata` in the following object is `a.vcf`:
        'outputs': {'merged': 'a.vcf'}

        :return: Filtered task outputs.
        """
        pass


class Metadata:
    """
    Implements utilities for traversing, processing, and
    querying the resulting metadata (in JSON) of running
    a workflow on Cromwell.
    """
    def __init__(self, filename):
        self.filename = filename

    @staticmethod
    def _get_output_label(parent_workflow, workflow, output_var, shard_index):
        """
        Composes a label for a task output.
        :return: Some examples of constructed labels are:
            - GATKSVPipelineSingleSample.Module00c.Module00c.PreprocessPESR.std_manta_vcf
            - Module00c.PreprocessPESR.PreprocessPESR.StandardizeVCFs.std_vcf.0
        """
        return \
            ((parent_workflow + ".") if parent_workflow else "") + \
            f"{workflow}.{output_var}" + \
            (("." + str(shard_index)) if shard_index != -1 else "")

    @staticmethod
    def _get_filtered_outputs(outputs):
        return outputs

    def _traverse_outputs(self, calls, parent_workflow="", deep=False):
        output_files = {}

        def update_output_files(outputs):
            if run["executionStatus"] == "Done" and len(outputs) > 0:
                output_files[self._get_output_label(
                    parent_workflow, workflow, out_label,
                    run["shardIndex"])] = outputs

        for workflow, runs in calls.items():
            for run in runs:
                if "outputs" in run:
                    for out_label, out_files in run["outputs"].items():
                        if not out_files:
                            continue
                        update_output_files(self._get_filtered_outputs(out_files))
                if deep and "subWorkflowMetadata" in run:
                    output_files.update(
                        self._traverse_outputs(
                            run["subWorkflowMetadata"]["calls"],
                            workflow, deep))
        return output_files

    def get_outputs(self, include_sub_workflows=False, filter_method=None):
        """
        Iterates through a given cromwell metadata file
        and filters the output files to be compared.

        :param include_sub_workflows: Boolean, if set to True,
        output files generated in sub-workflows will be traversed.

        :param filter_method: A method to override the default
        filter method. This method should be implement the
        ITaskOutputFilters interface. Every traversed output of tasks
        will be passed to this method, and this method's returned
        value will be aggregated and returned. For instance, see
        FilterBasedOnExtensions class for how the filter method can
        be used to extract files with certain extension from the
        metadata.

        :return: A dictionary with keys and values being a composite label
        for tasks outputs and the values of the task output, respectively.
        For instance (serialized to JSON and simplified for brevity):
        {
           "GATKSVPipelineSingleSample.FilterMelt.out":{
              "vcf.gz":["NA12878.melt.NA12878.vcf.gz"]
           }
        }
        """
        if filter_method:
            if not issubclass(type(filter_method.__self__),
                              ITaskOutputFilters):
                raise TypeError(f"The class {type(filter_method.__self__)} "
                                f"should implement the interface "
                                f"{ITaskOutputFilters}.")
            self._get_filtered_outputs = types.MethodType(filter_method, self)

        with open(self.filename, "r") as metadata_file:
            metadata = json.load(metadata_file)
            output_files = self._traverse_outputs(
                metadata["calls"],
                deep=include_sub_workflows)
        return output_files
