
# Bacterial Isolate Assembly (BIA) Pipeline for fully automated use.

For testing and manual use please have a look at https://github.com/seru71/bia_pipeline


## Pipeline

The pipeline consists of bcl2fastq, FastQC, bbmap, Trimmomatic, SPAdes, and QUAST


## Setup

Environment for the pipeline is set up using ansible playbook: https://github.com/seru71/ngs-compute-ansible

## Usage

* Starting the pipeline

	The NGS pipeline is started by a cron job when new data appears. 
	The exact commandline used is in the `monitor_new_datasets.sh`, which should be added to crontab.
	`monitor_new_datasets.sh` takes one argument - path to a folder where new runfolders are saved by the sequencer.

* Outputs

    The script will create RUN_ID folder in the scratch-root directory (given in pipeline config template). 
    The folder will contain: 
    - SAMPLE_ID/ - one dir per sample, named after samples listed in RUN_FOLDER/SampleSheet.csv and demultiplexed by bcl2fastq 
    - fastqs/    - raw FASTQ files created by bcl2fastq
    - drmaa/     - SLURM scripts created automatically (if you are using SLURM; for debugging purposes)
    - qc/        - qc output from FastQC and QUAST

    After finishing, the sample directories will contain:
    - FASTQ files at different processing stages
    - scaffolds in SAMPLE_ID/SAMPLE_ID_mra.fasta
    - contigs in SAMPLE_ID/SAMPLE_ID_mra_contigs.fasta
    - genome assembly files in SAMPLE_ID/SAMPLE_ID_mra (if not removed automatically)
    
    If archiving of fastqs/results is on (default; for details see pipeline.config.template), the fastqs/results & QC are copied to preconfigured archive location.
    The scratch directory can be cleared automatically by uncommenting `#@posttask(cleanup_files)` over `complete_run` definition in the pipeline, or by a cron job.
    
* Typical usage
	
	1. Play ansible playbook to set up the environment, paths, etc.
	The playbook adds `monitor_new_datasets.sh` to crontab with hourly executions.
	
	2. Sequence new samples storing runfolder with SampleSheet.csv in the location provided to `monitor_new_datasets.sh`
	
	3. Check the QC reports and final results in the archive. 
	If clearing the scratch dir hasn't been uncommented in the pipeline, intermediate results will be kept in the scratch dir.
	
	
    




