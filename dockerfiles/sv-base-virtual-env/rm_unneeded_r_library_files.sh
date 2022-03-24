#!/bin/bash
VIRTUAL_ENV_PATH=${1:-"/gatk-sv-env"}
R_LIB_PATH=$(find $VIRTUAL_ENV_PATH/R -name lib -type d | head -n1)

find $R_LIB_PATH \
    -type d \
    \( -name "help" -o -name "doc" -o -name "html" -o -name "htmlwidgets" -o -name "demo" -o -name "demodata" \
      -o -name "examples" -o -name "exampleData" -o -name "unitTests" -o -name "tests" -o -name "testdata" \
      -o -name "shiny" \) \
  | xargs rm -rf