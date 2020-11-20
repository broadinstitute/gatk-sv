#!/bin/bash

scripts/inputs/build_inputs.py input_templates inputs -a '{"ref_panel" : "ref_panel_1kg_v2"}'

scripts/inputs/build_inputs.py test_input_templates test_inputs_small -a '{"test_batch" : "test_batch_small"}'
scripts/inputs/build_inputs.py test_input_templates test_inputs_large -a '{"test_batch" : "test_batch_large"}'

scripts/inputs/build_inputs.py test_input_templates test_inputs_single_sample -a '{ "test_single_sample" : "test_single_sample_NA19240", "ref_panel" : "ref_panel_v1b" }'