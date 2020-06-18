#!/bin/bash
# Usage: watch_cromshell WORKFLOW_ID [sleep_time_sec]
# At regular polling interval, call cromshell-execution-status to get
# status of workflow / subworkflows. If any shards are failing their
# shard number will be listed. Terminates upon completion of workflow.
# -If cromshell is not located on your PATH, then define the ENV
#  variable CROMSHELL with the appropriate path.
# -Currently sub-sub-workflows are probably not handled correctly and
#  failing shards may incorrectly indicated in this scenario.

# having problems with this while debugging / improving. Eventually
# should reinstate:
#set -Eeuo pipefail

WORKFLOW_ID=$1
SLEEP_TIME=${2:-60}
CROMSHELL=${CROMSHELL:-"cromshell"}

WORKFLOW_RUNNING=true
FAILED_SHARDS=""
NUM_FAILED=0


function handle_status() {
    ID=$1
    DATE_STAMP=$(date)
    WF_STATUS_COUNT=$($CROMSHELL execution-status-count $ID 2> /dev/null)
    if [ -z "$WF_STATUS_COUNT" ] || [ "$WF_STATUS_COUNT" == "[]" ]; then
        # no status count is available. Possibly the job is just starting
        # up, or possibly there's an error
        WORKFLOW_STATUS="$($CROMSHELL status $ID 2>/dev/null | jq '."status"')"
        echo "WORKFLOW_STATUS=$WORKFLOW_STATUS"
        if [[ "$WORKFLOW_STATUS" =~ Submitted|Starting|Running ]]; then
            # if there are no status values, then the job status may
            # still be "Submitted". Check for that or "Running" (in case
            # status changes between the two checks. If it's either of
            # these, wait another sleep cycle to check status
            return 0
        else
            # the workflow is in an error state, but there is no data
            # about tasks. Probably an error with the WDL
            echo "Overall workflow failed with no execution status count"
            echo "Error messages:"
            METADATA=$($CROMSHELL metadata $WORKFLOW_ID 2>/dev/null)
            echo "$METADATA" | grep '"message":'
            WORKFLOW_RUNNING=false
            return 1
        fi
    fi
    NUM_TASKS=$(echo "$WF_STATUS_COUNT" | jq 'length')
    if [ $NUM_TASKS == 0 ]; then
        # haven't started yet, keep waiting
        return 0
    fi
    echo "$DATE_STAMP"
    NUM_TASKS_RUNNING=0
    for ((TASK_IND = 0; TASK_IND < NUM_TASKS; TASK_IND++)); do
        TASK_RUNNING=true
        TASK_STATUS_COUNT=$(echo "$WF_STATUS_COUNT" | jq ".[$TASK_IND]")
        check_status $ID "$TASK_STATUS_COUNT" $TASK_IND
        if $TASK_RUNNING; then
            ((NUM_TASKS_RUNNING++))
        fi
    done
    if [[ $NUM_TASKS_RUNNING == 0 ]]; then
        WORKFLOW_RUNNING=false
    else
        echo
    fi
}

function check_status() {
    ID=$1
    STATUS_COUNT="$2"
    TASK_NAME=$(echo "$STATUS_COUNT" | jq "keys | .[]")
    STATUS_VALUES=$(\
        echo "$STATUS_COUNT" \
        | grep -Ev "[]{}[]" | sed -e 's/[":,]//g' -e 's/^ *//' \
        | awk '{print $2 " " $1}' | paste -s -d, - | sed -e 's/,/, /g'\
    )
    NUM_WORKFLOWS=$(\
        echo "$STATUS_VALUES" \
            | awk ' BEGIN {
                    FS = "[, \t]+"
                } {
                for (i=1; i<NF; i+=2) {
                    STAT_VALUE = $(i+1)
                    if(STAT_VALUE != "RetryableFailure") {
                        TOT += $i
                    }
                }
            } END {
                print TOT
            }'\
    )
    echo "$TASK_NAME: $NUM_WORKFLOWS workflows: $STATUS_VALUES"
    
    if echo "$STATUS_VALUES" | grep -q "Failed"; then
        NUM_NEW_FAILED=$(echo "$STATUS_COUNT" | grep Failed | sed -E 's/[^0-9]*([0-9]*).*/\1/')
        if [ $NUM_NEW_FAILED != $NUM_FAILED ]; then
            # get metadata from this call
            METADATA=$(2>/dev/null $CROMSHELL metadata $WORKFLOW_ID)
            # get failed shard info by parsing metadata JSON:
            # 1) a) get status of each subworkflow as an array
            SUBWORKFLOW_STATUS=$( \
                echo "$METADATA" \
                | jq "
                    .calls.$TASK_NAME
                    | .[].subWorkflowMetadata.status
                  " \
            )
            if [ "$(echo "$SUBWORKFLOW_STATUS" | uniq)" == "null" ]; then
                # not calling subworkflows, calling tasks, so
                # 1) b) get status of each call as an array
                #    c) filtered out "Preempted" statuses
                CALL_STATUS=$( \
                    echo "$METADATA" \
                    | jq "
                        .calls.$TASK_NAME
                        | .[].backendStatus
                        | select(. != \"Preempted\")
                      " \
                )
            else
                CALL_STATUS=$SUBWORKFLOW_STATUS
            fi
            #   2) get the indices of "Failed" statuses as an array
            #   3) filter out calls that had no failed indices
            FAILED_SHARDS=$( \
                echo "$CALL_STATUS" \
                | jq -s '
                    indices("Failed")
                    | tostring
                  ' \
                | sed 's/[][]//g' \
            )
            NUM_FAILED=$NUM_NEW_FAILED
        fi
        if [[ ! -z "$FAILED_SHARDS" ]]; then
            echo "    Failed shards: $FAILED_SHARDS"
        fi
    fi
    
    if ! echo "$STATUS_VALUES" | grep -qE "Running|Starting|Submitted"; then
        # no more running, don't need to watch any longer
        TASK_RUNNING=false
    fi
}

echo "using workflow-id == $WORKFLOW_ID"
echo
handle_status $WORKFLOW_ID
while $WORKFLOW_RUNNING; do
    sleep $SLEEP_TIME
    handle_status $WORKFLOW_ID
done

# could give elapsed time with some effort
METADATA=${METADATA:-$(2>/dev/null $CROMSHELL metadata $WORKFLOW_ID)}
START_TIME=$(echo "$METADATA" | grep start | grep -v description \
             | cut -d'"' -f4 | sort | head -n1 \
             | sed -e 's/[A-z]$//' -e 's/[A-z]/ /g')
END_TIME=$(echo "$METADATA" | grep -E '("end"|"endTime")' \
          | cut -d'"' -f4 | sort | tail -n1 \
          | sed -e 's/[A-z]$//' -e 's/[A-z]/ /g')
echo "Start: $START_TIME"
echo "End:   $END_TIME"
