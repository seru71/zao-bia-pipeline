#!/bin/bash
#
# Script monitoring an input directory for new Illumina runs
#

# Directory to be monitored
DIR=$1

# Constants
REPO=https://github.com/seru71/zao-bia-pipeline.git
WORK_DIR=/ngs/prod/scratch
LOG_FILE=/tmp/monitor_new_dataset.log

# Check if new runfolder exists, and if the run is complete (marked with RTACommplete.txt)
# Looks for file RTAComplete with modification date <1h from now

new=`find $DIR/*/RTAComplete.txt -mmin -60`

echo `date -Iseconds`" - found: ${new}" >> ${LOG_FILE}
if [ ! -z "${new}" ]; then
	export RUN_FOLDER=`dirname ${new}`
	export RUN_ID=`basename ${RUN_FOLDER}`
	echo "Exporting RUN_FOLDER [${RUN_FOLDER}] and RUN_ID [${RUN_ID}]" >> ${LOG_FILE}
fi

# If the runfolder is complete, create pipeline config and start the pipeline
if [ ! -z $RUN_ID ]; then

	# Create pipeline config file
	git clone ${REPO} ${WORK_DIR}/${RUN_ID}/pipeline
#	sed 's///g' ${WORK_DIR}/${NEW_RUN}/pipeline/bia_pipeline.config.template > ${WORK_DIR}/${RUN_ID}/pipeline/bia_pipeline.config
	cp ${WORK_DIR}/${RUN_ID}/pipeline/bia_pipeline.config.template ${WORK_DIR}/${RUN_ID}/pipeline.config
	
	cmd="${WORK_DIR}/${RUN_ID}/pipeline/bia_pipeline.py \
		-s ${WORK_DIR}/${RUN_ID}/pipeline.config \
		--run_folder ${RUN_FOLDER} \
		-t complete_run \
		-j8 -vvv &> ${WORK_DIR}/${RUN_ID}/pipeline_`date -Iseconds`.err"

	echo -e "Starting the pipeline with command:\n${cmd}" >> ${LOG_FILE}

	# Start the pipeline
	${cmd}

fi
