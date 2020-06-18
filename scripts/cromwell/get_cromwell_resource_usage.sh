#!/bin/bash
# USAGE: get_cromwell_memory_usage.sh WORKFLOW_ID
#              or
#        get_cromwell_memory_usage.sh GCS_PATH_TO_WORKFLOW_FOLDER
# Displays #shards x #fields table of memory usage info with one header
# line.
# -If you pass a workflow id, cromshell will be used to find the
#  appropriate workflow folder in google cloud storage.
# -If cromshell is not located on your PATH, then define the ENV
#  variable CROMSHELL with the appropriate path.
# -This script works by finding all the logs, sorting them into sensible
#  order, and chopping up their paths to make a description column. If
#  any jobs have not completed they will simply be omitted. The script
#  makes no attempt to figure out what tasks *should* have run. However,
#  the description field should make any such omissions discoverable.
# -Since there is significant start-up time to running gsutil functions,
#  running the inner loop of this script in parallel results in
#  signficant speed-up. Installing parallel (available on osx and linux)
#  triggers this automatically.

set -Eeu -o pipefail

CROMSHELL=${CROMSHELL:-"cromshell"}

WORKFLOW_INFO="$1"
if [[ $WORKFLOW_INFO == "gs://"* ]]; then
    WORKFLOW_DIR="$WORKFLOW_INFO"
else
    WORKFLOW_ID="$WORKFLOW_INFO"
    # get the metadata for this workflow id
    METADATA=$(2>/dev/null "$CROMSHELL" metadata $WORKFLOW_ID)
    # get the appropriate path in google cloud for the workflow dir
    # from the metadata
    # a) find lines in the metadata that include a 
    WORKFLOW_DIR=$( \
        echo "$METADATA" \
        | grep -o "gs://[^\S]*$WORKFLOW_ID" \
        | tail -n1 \
    )
fi
1>&2 echo "WORKFLOW_DIR=$WORKFLOW_DIR"


function get_monitor_logs() {
    TOP_DIR=$1
    gsutil -m ls "$TOP_DIR/**monitoring.log" 2>/dev/null || echo ""
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
              -v SHARD_FORMAT="$SHARD_FORMAT" ' {
            if($1 == "shard") {
                SHARD_KEY=sprintf("%s/" SHARD_FORMAT, SHARD_KEY, $2)
            } else if($1 == "call") {
                CALL_KEY=sprintf("%s/%s", CALL_KEY, $2)
            } else if($1 == "attempt") {
                ATTEMPT_NUMBER=$2
            }
          } END {
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
        | while read LOG_FILE; do
            SORT_KEY=$(get_task_sort_key "$LOG_FILE" "$TOP_DIR" "$NUM_LOGS")
            printf "%s\t%s\n" "$LOG_FILE" "$SORT_KEY"
          done \
        | sort -t $'\t' -k3,3 -k2,2rn \
        | uniq -f2 \
        | cut -d$'\t' -f1
}


function date_stamp_to_seconds() {
    if [ "$(uname)" == "Darwin" ]; then
        date -j -u -f "%a %b %d %T %Z %Y" "$1" "+%s"
    else
        date -d "$1" "+%s"
    fi
}
export -f date_stamp_to_seconds


function get_task_peak_resource_usage() {
    LOG_FILE=$1
    
    gsutil cat "$LOG_FILE" \
        | awk -v OFS='\t' '
            /^\[.*\]$/ {
                DATE_STR=substr($0,2,length($0)-2)
                if(MIN_TIME == "") {
                    MIN_TIME=DATE_STR
                }
                else {
                    MAX_TIME=DATE_STR
                }
            }
            $1 == "*" {
                if($2 == "Memory") {
                    if($4 > PEAK_MEM) {
                        PEAK_MEM = $4
                    }
                } else if($2 == "CPU") {
                    if($4 > PEAK_CPU) {
                        PEAK_CPU = $4
                    }
                } else if($2 == "Disk") {
                    DISK=$4
                    LEN=length(DISK)
                    UNIT=substr(DISK, LEN)
                    if(UNIT ~ /[A-Z]/) {
                        if(UNIT == "T") {
                            SCALE=2^10
                        } else if(UNIT == "G") {
                            SCALE=1
                        } else if(UNIT == "M") {
                            SCALE=2^-10
                        } else if(UNIT == "K") {
                            SCALE=2^-20
                        } else if(UNIT == "B") {
                            SCALE=2^-30
                        } else {
                            SCALE=1
                        }
                        DISK_VAL=substr(DISK, 0, LEN-1) * SCALE
                    } else {
                        DISK_VAL=DISK
                    }
                    if(DISK_VAL > PEAK_DISK) {
                        PEAK_DISK=DISK_VAL
                    }
                } else if($2 == "Read/Write") {
                    if($4 > PEAK_READ || PEAK_READ == "") {
                        PEAK_READ=$4
                    }
                    if($6 > PEAK_WRITE || PEAK_WRITE == "") {
                        PEAK_WRITE=$6
                    }
                }
            } END {
                if(PEAK_MEM == "") {
                    PEAK_MEM="nan"
                }
                if(PEAK_DISK == "") {
                    PEAK_DISK="nan"
                }
                if(PEAK_CPU == "") {
                    PEAK_CPU="nan"
                }
                if(PEAK_READ == "" || PEAK_READ == "N/A") {
                    PEAK_READ="nan"
                }
                if(PEAK_WRITE == "" || PEAK_WRITE == "N/A") {
                    PEAK_WRITE="nan"
                }
                if(MAX_TIME == "") {
                    MIN_TIME = "Thu Jan 1 00:00:00 UTC 1970"
                    MAX_TIME = MIN_TIME
                }
                print PEAK_MEM, PEAK_DISK, PEAK_CPU, PEAK_READ, PEAK_WRITE
                print MIN_TIME
                print MAX_TIME
            }'
}
export -f get_task_peak_resource_usage


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
    DESCRIPTION=$(get_task_description "$LOG_FILE" "$TOP_DIR")
    RESOURCE_USAGE=$(get_task_peak_resource_usage "$LOG_FILE")
    get_task_peak_resource_usage "$LOG_FILE" | {
        read RESOURCE_USAGE
        read MIN_TIME
        read MAX_TIME
        MIN_SECONDS=$(date_stamp_to_seconds "$MIN_TIME")
        MAX_SECONDS=$(date_stamp_to_seconds "$MAX_TIME")
        awk -v DESCRIPTION="$DESCRIPTION" \
            -v RESOURCE_USAGE="$RESOURCE_USAGE" \
            -v MIN_SECONDS=$MIN_SECONDS \
            -v MAX_SECONDS=$MAX_SECONDS \
            'END {
                printf "%s\t%.3f\t%s\n", RESOURCE_USAGE, (MAX_SECONDS-MIN_SECONDS)/60/60, DESCRIPTION 
            }' /dev/null
    }
}
export -f get_task_columns


function get_workflow_peak_resource_usage() {
    export TOP_DIR=$1
    echo -e "mem_GiB\tdisk_GiB\tcpu_%\tread_MiB/s\twrite_MiB/s\truntime_Hours\ttask_description"
    LOGS=$(get_monitor_logs "$TOP_DIR" | sort_monitor_logs "$TOP_DIR")
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
            BAR="--bar"
        fi
        echo "$LOGS" | parallel ${BAR} --env TOP_DIR -k "get_task_columns {} $TOP_DIR"
    else
        1>&2 echo "Consider installing 'parallel', it will give significant speed-up"
        echo "$LOGS" | while read WORKFLOW_LOG; do
            get_task_columns "$WORKFLOW_LOG" "$TOP_DIR"
        done
    fi
}

get_workflow_peak_resource_usage "$WORKFLOW_DIR"
