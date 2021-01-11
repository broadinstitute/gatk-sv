#!/bin/bash

function show_help() {
    cat <<-END
USAGE: get_cromwell_memory_usage2.sh [OPTIONS] WORKFLOW_INFO
Displays #tasks x #fields table of resource usage info with two
header lines, and additional column of task descriptions.
    WORKFLOW_INFO: specifies workflow identity. Can be
                   a) a cromwell workflow ID
                   b) a path to workflow output in google cloud (starting with gs://)
                   c) a local path to workflow output
    OPTIONS:
        -r --raw-output
            If set, output data as tab-separated table. Otherwise pass data through column -t for
            easier reading.
        -u --no-units
            If set, don't show second header line with units.
        -o --output-file OUTPUT_FILE_NAME
            If specified, write output to file instead of stdout. Note that all non-tabular output is
            sent to stderr, so this is just syntactic sugar for file redirection.

-This script works by finding all the logs, sorting them into sensible
 order, and chopping up their paths to make a description column. If
 any jobs have not completed they will simply be omitted. The script
 makes no attempt to figure out what tasks *should* have run. However,
 the description field should make any such omissions discoverable.
-If you run on a local path, the log file names must still be
 "monitoring.log", and the local folder structure must be the same as
 in the original cloud bucket (other non-log files are not required
 though)
-If you pass a workflow id, cromshell will be used to find the
 appropriate workflow folder in google cloud storage.
-If cromshell is not located on your PATH, then define the ENV
 variable CROMSHELL with the appropriate path.
-Since there is significant start-up time to running gsutil functions,
 running the inner loop of this script in parallel results in
 signficant speed-up. Installing gnu parallel (available on osx and
 linux) triggers this automatically.
END
}

if [[ $# == 0 ]]; then
    show_help
    exit 0
fi
RAW_OUTPUT=false
SHOW_UNITS=true
OUTPUT_FILE="/dev/stdout"
WORKFLOW_INFO=""
for ((i=1; i<=$#; ++i)); do
    if [[ ${!i} =~ ^-+(h|help) ]]; then
        show_help
        exit 0
    elif [[ ${!i} =~ ^-+(r|raw-output) ]]; then
        RAW_OUTPUT=true
    elif [[ ${!i} =~ ^-+(u|no-units) ]]; then
        SHOW_UNITS=false
    elif [[ ${!i} =~ ^-+(o|output-file) ]]; then
        ((++i))
        OUTPUT_FILE="${!i}"
    elif [[ ${!i} =~ ^-.* ]]; then
        1>&2 echo "Unknown option ${!i}"
        show_help
        exit 1
    elif [[ -z "$WORKFLOW_INFO" ]]; then
        WORKFLOW_INFO=${!i}
    else
        1>&2 echo "Too many arguments"
        exit 1
    fi
done
export SHOW_UNITS

set -Eeu -o pipefail

CROMSHELL=${CROMSHELL:-"cromshell"}

REMOTE=true
if [[ -z "$WORKFLOW_INFO" ]]; then
    1>&2 echo "No WORKFLOW_INFO provided"
    show_help
    exit 1
elif [[ $WORKFLOW_INFO == "gs://"* ]]; then
    # workflow info is a cloud file, strip trailing slashes (must use sed because OSX uses old bash)
    WORKFLOW_DIR=$(echo "$WORKFLOW_INFO" | sed 's,/*$,,g')
elif [[ -d "$WORKFLOW_INFO" ]]; then
    # workflow info is a local file, strip trailing slashes (must use sed because OSX uses old bash)
    WORKFLOW_DIR=$(echo "$WORKFLOW_INFO" | sed 's,/*$,,g')
    REMOTE=false
else
    WORKFLOW_ID="$WORKFLOW_INFO"
    # get the metadata for this workflow id
    if ! METADATA=$(2>/dev/null $CROMSHELL -t 60 slim-metadata $WORKFLOW_ID); then
        1>&2 echo "Unable to obtain workflow $WORKFLOW_ID metadata from cromshell, try supplying GCS_PATH_TO_WORKFLOW_FOLDER"
        exit 1
    fi
    # get the appropriate path in google cloud for the workflow dir
    # from the metadata
    # a) find lines in the metadata that include a 
    WORKFLOW_DIR=$( \
        echo "$METADATA" \
        | grep -Eo "gs://[^[:space:]]*$WORKFLOW_ID" \
        | tail -n1 \
    )
fi
1>&2 echo "WORKFLOW_DIR=$WORKFLOW_DIR"


function get_monitor_logs() {
    TOP_DIR=$1
    REMOTE=$2
    if $REMOTE; then
      gsutil -m ls "$TOP_DIR/**monitoring.log" 2>/dev/null || echo ""
    else
      find "$TOP_DIR" -name "monitoring.log" 2>/dev/null || echo ""
    fi
}


# ingest LOG_FILE, TOP_DIR, NUM_SHARDS
# print out sorting key of form ATTEMPT_KEY -tab- MAIN_KEY
#    where ATTEMPT_KEY
#       -lists the attempt number of the executing tasks OR
#       -if there is no attempt in the path, calls it attempt 0
#       -digits are padded the same as shards
#    where MAIN_KEY
#       -preserves info about calls and shard numbers, each separated
#        by '/'
#       -shard numbers having enough 0-padded digts to all be of
#        the same length
function get_task_sort_key() {
    LOG_FILE=$1
    TOP_DIR=$2
    N_START=$((1 + $(echo "$TOP_DIR" | tr / '\n' | wc -l)))
    NUM_SHARDS=$(($3))
    MAX_SHARD_DIGITS=${#NUM_SHARDS}
    SHARD_FORMAT="%0${MAX_SHARD_DIGITS}d"
    # keep info about task calls, shards, and attempts below top-dir
    # if there is no preemption folder in the path, call it attempt 0
    echo "$LOG_FILE" \
        | tr / '\n' \
        | tail -n+$N_START \
        | awk -v FS='-' \
              -v SHARD_FORMAT="$SHARD_FORMAT" '
            {
                if($1 == "shard") {
                    SHARD_KEY=sprintf("%s/" SHARD_FORMAT, SHARD_KEY, $2)
                } else if($1 == "call") {
                    CALL_KEY=sprintf("%s/%s", CALL_KEY, $2)
                } else if($1 == "attempt") {
                    ATTEMPT_NUMBER=$2
                }
            }
            END {
              printf SHARD_FORMAT "\t%s/%s", ATTEMPT_NUMBER, CALL_KEY, SHARD_KEY
            }'
}


function sort_monitor_logs() {
    TOP_DIR="$1"
    LOGS_LIST=$(cat)
    NUM_LOGS=$(($(echo "$LOGS_LIST" | wc -l)))
    # The older bash on OSX does not have associative arrays, so to
    # sort file names according to a key, we join the key and the file
    # name into one string with tab delimiters (okay because these are
    # cloud paths produced by cromwell and have no tabs). Then sort by
    # the key, and ultimately cut away the key. There is one extra
    # complication that there may be multiple "attempts" at each task,
    # and we only want to keep the final (presumably successful)
    # attempt.
    #
    # 1. for each log file
    #  a) get a sort key of form: ATTEMPT_KEY tab MAIN_KEY
    #  b) print line of form: LOG_FILE tab MAIN_KEY tab ATTEMPT_KEY
    # 2. sort lines by increasing MAIN_KEY, and secondarily by
    #      decreasing (numeric) ATTEMPT_KEY
    # 3. keep first unique instance of MAIN_KEY (i.e. the last attempt)
    # 4. print out the log file (the first field) in sorted order
    echo "$LOGS_LIST" \
        | while read -r LOG_FILE; do
            SORT_KEY=$(get_task_sort_key "$LOG_FILE" "$TOP_DIR" "$NUM_LOGS")
            printf "%s\t%s\n" "$LOG_FILE" "$SORT_KEY"
          done \
        | sort -t $'\t' -k3,3 -k2,2rn \
        | uniq -f2 \
        | cut -d$'\t' -f1
}


# Scan LOG_FILE, extract header, and print maximum of each column.
# If a column is missing data, print "nan"
function get_task_peak_resource_usage() {
    LOG_FILE=$1
    if $REMOTE; then
        gsutil cat "$LOG_FILE"
    else
        cat "$LOG_FILE"
    fi \
        | awk '
            function handle_nan(num) {
               return num < 0 ? "nan" : num
            }
            BEGIN {
                NEED_HEADER=2
            }
            NEED_HEADER == 0 {
                split($0, WORDS, /\t/)
                PEAK_VALUE[1] = WORDS[1]
                for(i=2; i<=length(WORDS); ++i) {
                    WORD=WORDS[i]
                    if(length(WORD) > 0 && WORD > PEAK_VALUE[i]) {
                        PEAK_VALUE[i] = WORD
                    }
                }
            }
            NEED_HEADER>0 {
                if(NEED_HEADER==2) {
                    if($1" "$2 == "Num processors:") {
                        TOT["CPU"] = $3
                        TOT_NAME["CPU"] = "nCPU"
                        TOT_UNIT["CPU"] = "#"
                    }
                    else if($1" "$2 == "Total Memory:") {
                        TOT["Mem"] = $3
                        TOT_NAME["Mem"] = "TotMem"
                        TOT_UNIT["Mem"] = $4
                    }
                    else if($1" "$2 == "Total Disk") {
                        TOT["Disk"] = $4
                        TOT_NAME["Disk"] = "TotDisk"
                        TOT_UNIT["Disk"] = $5
                    }
                    else if($1 == "ElapsedTime") {
                        # this is the first header line
                        # for summary purposes, augment instantaneous
                        # usage with total VM stats
                        PEAK_VALUE[1] = "00:00:00"
                        HEADER_NAME[1] = $1
                        printf "%s", $1
                        for(i=2; i<=NF; ++i) {
                            PEAK_VALUE[i] = -1.0
                            HEADER_NAME[i] = $i
                            if($i in TOT_NAME) {
                              printf "\t%s", TOT_NAME[$i]
                              delete TOT_NAME[$i]
                            }
                            printf "\t%s", $i
                        }
                        printf "\n"
                        --NEED_HEADER
                    }
                } else {
                    printf "%s", $1
                    for(i=2; i<=NF; ++i) {
                        NAME_i = HEADER_NAME[i]
                        if(NAME_i in TOT_UNIT) {
                            printf "\t%s", TOT_UNIT[NAME_i]
                            delete TOT_UNIT[NAME_i]
                        }
                        printf "\t%s", $i
                    }
                    printf "\n"
                    --NEED_HEADER
                }
            }
            END {
                printf "%s", handle_nan(PEAK_VALUE[1])
                for(i=2; i<=length(PEAK_VALUE); ++i) {
                    NAME_i = HEADER_NAME[i]
                    if(NAME_i in TOT) {
                        printf "\t%s", handle_nan(TOT[NAME_i])
                        delete TOT[NAME_i]
                    }
                    printf "\t%s", handle_nan(PEAK_VALUE[i])
                }
                printf "\n"
            }
        '
}
export -f get_task_peak_resource_usage


# Condense directory structure of full path to LOG_FILE into a succinct
# description of the task. Ignore components above TOP_DIR, as they are
# common to all the log files that are being requested.
function get_task_description() {
    LOG_FILE=$1
    if [ $# -ge 2 ]; then
        TOP_DIR=$2
        N_START=$((1 + $(echo "$TOP_DIR" | tr / '\n' | wc -l)))
    else
        N_START=1
    fi
    # keep info about task calls and shards below top-dir
    echo "$LOG_FILE" \
        | tr / '\n' \
        | tail -n+$N_START \
        | grep -E "^(call-|shard-|attempt-)" \
        | tr '\n' / \
        | sed -e 's/call-//g' -e 's,/$,,'
}
export -f get_task_description

function get_task_columns() {
    LOG_FILE="$1"
    TOP_DIR="$2"
    RESOURCE_USAGE=$(get_task_peak_resource_usage "$LOG_FILE")
    if [[ -n "$RESOURCE_USAGE" ]]; then
        DESCRIPTION=$(get_task_description "$LOG_FILE" "$TOP_DIR")
    
        # due to OSX having an ancient version of bash, this produces syntax errors:
        # paste <(echo "$RESOURCE_USAGE" | head -n2) <(echo "task")
        printf "%s\ttask\n" "$(echo "$RESOURCE_USAGE" | head -n1)"
        if $SHOW_UNITS; then
          echo "$RESOURCE_USAGE" | tail -n2 | head -n1
        fi
        printf "%s\t%s\n" "$(echo "$RESOURCE_USAGE" | tail -n1)" "$DESCRIPTION"
    fi
}
export -f get_task_columns


function get_workflow_peak_resource_usage() {
    export TOP_DIR=$1
    export REMOTE=$2
    LOGS=$(get_monitor_logs "$TOP_DIR" $REMOTE | sort_monitor_logs "$TOP_DIR")
    if [ -z "$LOGS" ]; then
        1>&2 echo "No logs found in $TOP_DIR"
        exit 0
    fi
    if command -v parallel > /dev/null; then
        # parallel command is installed, use it, much faster!
        if [ -t 1 ]; then
            # stdout is a terminal, not being redirected, don't use bar
            BAR=""
        else
            # being redirected, show progress via bar to stderr
            #BAR="--bar"
            # NOTE: keeping above line for now to see if I can find a way
            # to make it work, but it looks like --bar may be incompatible
            # with filtering out potentially empty results from logs that
            # were truncated before the header line
            BAR=""
        fi
        
        echo "$LOGS" \
            | parallel ${BAR} --env TOP_DIR -k --colsep $'\t' "get_task_columns {1} $TOP_DIR"
    else
        1>&2 echo "Consider installing 'parallel', it will give significant speed-up"
        echo "$LOGS" | while read -r WORKFLOW_LOG; do
            get_task_columns "$WORKFLOW_LOG" "$TOP_DIR"
        done
    fi \
        | if $SHOW_UNITS; then
            awk 'FNR < 3 || FNR%3 == 0 { print $0 }'
          else
            awk 'FNR < 2 || FNR%2 == 0 { print $0 }'
          fi
}

if $RAW_OUTPUT; then
    get_workflow_peak_resource_usage "$WORKFLOW_DIR" $REMOTE
else
    get_workflow_peak_resource_usage "$WORKFLOW_DIR" $REMOTE | column -t
fi > "$OUTPUT_FILE"
