
# Bacterial Isolate Assembly (BIA) Pipeline 



## Pipeline

The pipeline consists of bcl2fastq, Trimmomatic, FastQC, and SPAdes



## Setup

### Pipeline

1. Install Python if you don't have it from before, and a python package called `Ruffus` (http://www.ruffus.org.uk/). 
Running jobs on a cluster (PBS, Slurm, etc) requires `drmaa` package. 
You might also need following packages: optparse, logging, shutil

2. Clone the pipeline repository:
`git clone https://github.com/seru71/bia_pipeline.git <PIPELINE_HOME>`

3. (optional) Change directory to newly created pipeline dir, and checkout desired version (master by default)
```
cd <PIPELINE_HOME>
git checkout v0.1
```

The pipeline is ready now, but you will need all of its components to perform the analysis.

### Pipeline components

Install following tools:
1. FastQC (https://www.bioinformatics.babraham.ac.uk/projects/fastqc)
2. BBMap (https://sourceforge.net/projects/bbmap)
2. Trimmomatic (http://www.usadellab.org/cms/?page=trimmomatic)
3. SPAdes (http://bioinf.spbau.ru/spades)
4. QUAST (http://bioinf.spbau.ru/quast)

And if your input is not yet converted to FASTQ, install also:  
6. bcl2fastq ()


### Reference data

Download reference genome of the bacteria you are planning to analyze.

## Usage

The NGS pipeline is run using `bia_pipeline.py` script:

* Running the script

    You can run the script using `python <PIPELINE_HOME>/bia_pipeline.py`.
    A list of possible options will be presented. 
    Typically you will need to specify config/settings file (`-s`) and target task to be executed (`-t`). 
    If bcl2fastq conversion is to be run, you will need to specify run folder (`--run_folder RUN_FOLDER`) which indicates location of the run folder created by Illumina sequencer.
    Other useful arguments are dry-run (`-n`), verbosity level (`-v`), and number of concurrent jobs (`-j`).
    
    Important part of the pipeline is the config file which contains paths to tools, reference genome, and docker settings.
    See an exemplary file for all required options and documentation in `<PIPELINE_HOME>/bia_pipeline.config`
    If the settings file is not given as argument ( `--settings` / `-s` ), it is expected in the `<RUN_FOLDER>/bia_pipeline.config`
  
    Currently the pipeline provides two ways of preparing PE reads for the assembly.
    The default is to first merge overlapping pairs using bbmerge, trim both merged and unmerged reads, and use in the assembly merged reads, trimmed R1 and R2, as well as trimmed unpaired R1.
    This approach generates a merged-read assembly (mr_assembly folder) for every sample.
    To use the default read-merging approach, use `-t complete_run`.
    
    Another option is using trimmed PE reads without merging (trimmed-reads assembly; tr_assembly folder).
    Then bbmerge is not used, and SPAdes is provided with trimmed R1 and R2 FASTQs, and trimmed unpaired R1.
    To use the alternative approach, use `-t assemble_trimmed` or `-t qc_tr_assemblies`.
    
    The two approaches can be run side by side and the results compared.
    To compare across both types of assemblies run QUAST on `<scratch-root>/<RUN_ID>/*/[mt]r_assembly/contigs.fasta`
  
    If you want to follow progress of the pipeline script, use the verbose option (`-vvvvv`).
    In order to use multithreading, use the `-j` option (`-j 12`).

* Outputs

    The script will create RUN_ID folder in the scratch-root directory (given in settings). 
    Inside there will be several directories: 
    - SAMPLE_ID/ - one dir per sample, named after samples found in the RUN_FOLDER 
    - fastqs/    - fastq files
    - drmaa/     - SLURM scripts created automatically (if you are using SLURM; for debugging purposes)
    - qc/        - qc output

    After finishing, the sample directories will contain:
    - FASTQ files at different processing stages
    - assembled genome
    - FastQC reports on raw and processed FASTQs
    - QUAST report on the assemblies
   

* Typical usage

    For running the assembly using 12 concurrent threads:

	`bia_pipeline.py --run_folder /incoming/RUN_XXXX_YYYY_ZZZZZ \
    
						    --settings /my_project/ecoli.config \
                            
							--target complete_run \
                            
							-vvv -j 12 \
                            
							--log_file /my_project/pipeline_run_XXX.log`





