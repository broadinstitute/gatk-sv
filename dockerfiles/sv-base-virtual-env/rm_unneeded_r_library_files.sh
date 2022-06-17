#!/bin/bash
R_INSTALL_PATH=${1:-/opt/R}
R_LIB_PATH=$(find $R_INSTALL_PATH -name lib -type d -print -quit)

find $R_LIB_PATH \
    -type d \
    \( -name "help" -o -name "doc" -o -name "html" -o -name "htmlwidgets" -o -name "demo" -o -name "demodata" \
      -o -name "examples" -o -name "exampleData" -o -name "unitTests" -o -name "tests" -o -name "testdata" \
      -o -name "shiny" \) \
  | xargs rm -rf