#!/bin/bash
sed -e 's/\s/ /g' -e 's/  */ /g' $@
