#!/bin/bash
#
# Script monitoring an input directory for new Illumina runs
#

# Directory to be monitored
DIR=$1

# Constants
REPO=https://github.com/seru71/zao-bia-pipeline.git


# Check if new folders exist, and if the run is complete
# Looks for file RTAComplete with mod date <1h from now

touch /tmp/monitor_new_dataset.stamp
new=`find $DIR/RTAComplete.txt -mmin -60`

if [ ! -z "$new" ]; then
	export RUN_FOLDER=`dirname $new`
	export RUN_ID=`basename $RUN_FOLDER`
fi

# If the runfolder is complete, create pipeline config and start the pipeline
if [ ! -z $RUN_ID ]; then

	# Create pipeline config file
	git clone ${REPO} ${WORK_DIR}/${RUN_ID}/pipeline
	sed 's///g' ${WORK_DIR}/${NEW_RUN}/pipeline/bia_pipeline.config.template > ${WORK_DIR}/${RUN_ID}/pipeline/bia_pipeline.config

	# Start the pipeline
	{WORK_DIR}/${RUN_ID}/pipeline/bia_pipeline.py \
		-s ${WORK_DIR}/${RUN_ID}/pipeline/bia_pipeline.config \
		--run_folder ${NEW_RUN} \
		-t complete_run \
		-j8 -vvv &> {WORK_DIR}/${RUN_ID}/pipeline_`date -Iseconds`.err

fi
