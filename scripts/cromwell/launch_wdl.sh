#!/bin/bash
set -eu -o pipefail

WDL=$(realpath $1)
CONFIG_FILE=${2:-"$HOME/code/cromwell/cromwell_workflow_options.json"}
VALIDATE=${VALIDATE:-false}


SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
GATK_SV_ROOT=$SCRIPT_DIR
while [ $(basename "$GATK_SV_ROOT") != gatk-sv ]; do
    GATK_SV_ROOT=$(dirname "$GATK_SV_ROOT")
done

WDL_FILENAME=$(basename "$WDL")
WDL_NAME=${WDL_FILENAME%.*}

CLOUD_ENV="$GATK_SV_ROOT/inputs/values/google_cloud.my_project.json"
echo "CLOUD_ENV=$CLOUD_ENV"
cat << EOF > "$CLOUD_ENV"
{
  "google_project_id": "broad-dsde-methods",
  "terra_billing_project_id": "broad-dsde-methods"
}
EOF


RUN_DIR="$GATK_SV_ROOT/runs/$WDL_NAME"
DEPS_ZIP="$RUN_DIR/deps.zip"
rm -rf "$RUN_DIR"
mkdir -p "$RUN_DIR"
cd "$(dirname $WDL)"
zip "$DEPS_ZIP" *.wdl &> /dev/null
cd "$GATK_SV_ROOT"
"$GATK_SV_ROOT/scripts/inputs/build_default_inputs.sh" \
  -d "$GATK_SV_ROOT" \
  -c google_cloud.my_project \
  > /dev/null

rm -f $CLOUD_ENV

echo "Available input jsons:"
printf "%d\t%s\n" 0 "none (skip cromwell submit)"
n=1
JSON_ARRAY=()
while read INPUT_JSON; do
  # printf "%d\t%s\n" $n $(basename "$INPUT_JSON")
  printf "%d\t%s\n" $n "$INPUT_JSON"
  JSON_ARRAY["$n"]="$INPUT_JSON"
  n=$((n+1))
done < <(find "$GATK_SV_ROOT/inputs/build" -name "$WDL_NAME"*.json | sort)

read -p "Select input json number: " INPUT_JSON_NUMBER
if [[ $INPUT_JSON_NUMBER == 0 ]]; then
  exit 0
fi
INPUT_JSON=${JSON_ARRAY[$INPUT_JSON_NUMBER]}
cp "$INPUT_JSON" "$RUN_DIR"

if $VALIDATE; then
  womtool validate  $WDL -i "$INPUT_JSON"
else
  cromshell submit "$WDL" "$INPUT_JSON" "$CONFIG_FILE" "$DEPS_ZIP"
fi

