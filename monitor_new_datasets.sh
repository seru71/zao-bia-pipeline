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


# Check directory to be monitored
if [ -z $1 -o ! -e $1 ]; then
    log "Missing or incorrect path to the directory to be monitored."
    exit(1)
fi

export DIR=$1


# Check if new runfolder exists, and if the run is complete (marked with RTAComplete.txt)
# Looks for file RTAComplete with modification date <1h from now
new=`find $DIR/*/RTAComplete.txt -mmin -60`

if [ -z "${new}" ]; then
	log "No new runs found in $DIR"
    exit(0)
fi

log "Found: ${new}"

# Configure workspace for the pipeline
log "Configuring workspace for RUN_FOLDER [${RUN_FOLDER}] and RUN_ID [${RUN_ID}]"
export RUN_FOLDER=`dirname ${new}`
export RUN_ID=`basename ${RUN_FOLDER}`

git clone ${REPO} ${WORK_DIR}/${RUN_ID}/pipeline
cp ${WORK_DIR}/${RUN_ID}/pipeline/bia_pipeline.config.template ${WORK_DIR}/${RUN_ID}/pipeline.config
	
cmd="${WORK_DIR}/${RUN_ID}/pipeline/bia_pipeline.py \
-s ${WORK_DIR}/${RUN_ID}/pipeline.config \
--run_folder ${RUN_FOLDER} \
-t complete_run \
-j8 -vvv &> ${WORK_DIR}/${RUN_ID}/pipeline_`date -Iseconds`.err"


# Start the pipeline
log "Starting the pipeline with command: ${cmd}"
${cmd}
log "Pipeline finished with exit code: $?"
