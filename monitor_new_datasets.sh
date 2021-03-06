#!/bin/bash
#
# Script monitoring an input directory for new Illumina runs
#

# Constants
REPO=https://github.com/seru71/zao-bia-pipeline.git
WORK_DIR=/ngs/prod/scratch
LOG_FILE=/tmp/monitor_new_dataset.log

function log {
    echo -e `date -Iseconds`" - $1"
}

function run_as_ngs {
    echo "$*" | sudo -u ngs -i bash
}

function update_myself {
    # update myself
    cwd=`dirname "$(readlink -f "$0")"`
    cd $cwd
    # clean local changes
    git checkout -f master
    res=`git pull`
    if [ "$res" != "Already up-to-date." ]; then
		# if updated, run itself and exit with the same exit code
		./monitor_new_datasets.sh "$*"
		exit $?
	fi
}

update_myself "$*"

# Check directory to be monitored
if [ -z "$1" ] || [ ! -d $1 ]; then
    log "Missing or incorrect path to the directory to be monitored: [$1]"
    exit 1
fi

export DIR=$1


# Check if new runfolder exists, and if the run is complete (marked with RTAComplete.txt)
# Looks for file RTAComplete with modification date <1h from now
new=`find $DIR/*/RTAComplete.txt -mmin -60`

if [ -z "${new}" ]; then
    log "No new runs found in ${DIR}"
    exit 0
fi

log "Found: ${new}"

# Configure workspace for the pipeline
export RUN_FOLDER=`dirname ${new}`
export RUN_ID=`basename ${RUN_FOLDER}`
log "Configuring workspace for RUN_FOLDER [${RUN_FOLDER}] and RUN_ID [${RUN_ID}]"

if [ -e ${WORK_DIR}/${RUN_ID}/pipeline ]; then
    log "Found existing pipeline directory in [${WORK_DIR}/${RUN_ID}/pipeline]. Exiting."
    exit 2
fi

run_as_ngs git clone ${REPO} ${WORK_DIR}/${RUN_ID}/pipeline
run_as_ngs cp ${WORK_DIR}/${RUN_ID}/pipeline/bia_pipeline.config.template ${WORK_DIR}/${RUN_ID}/pipeline.config
	
timestamp=`date -Iseconds`
cmd="nice ${WORK_DIR}/${RUN_ID}/pipeline/bia_pipeline.py \
-s ${WORK_DIR}/${RUN_ID}/pipeline.config \
--run_folder ${RUN_FOLDER} \
-t complete_run \
-j8 -vvv -L ${WORK_DIR}/${RUN_ID}/pipeline_${timestamp}.log &> ${WORK_DIR}/${RUN_ID}/pipeline_${timestamp}.err"


# Start the pipeline
log "Starting the pipeline for $RUN_ID with command: ${cmd}"
run_as_ngs "${cmd}"
log "Pipeline for run ${RUN_ID} finished with exit code: $?"

exit 0

