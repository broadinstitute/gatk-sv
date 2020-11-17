#!/bin/bash

scripts/inputs/build_inputs.py input_templates inputs -a '{"ref_panel" : "ref_panel_1kg"}'

scripts/inputs/build_inputs.py test_input_templates test_inputs_small -a '{"test_batch" : "test_batch_small", "ref_panel" : "ref_panel_v1b"}'