#!/usr/bin/env python
"""

    bia_pipeline.py
                        --run_folder PATH
                        [--settings PATH] 
                        [--log_file PATH]
                        [--verbose]
                        [--target_tasks]  (by default the last task in the pipeline)
                        [--jobs N]        (by default 1)
                        [--just_print]
                        [--flowchart]
                        [--key_legend_in_graph]
                        [--forced_tasks]
                        [--run_on_bcl_tile TILE_REGEX]

"""
import sys
import os
import glob


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   options


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888


if __name__ == '__main__':
    from optparse import OptionParser
    import StringIO

    parser = OptionParser(version="%prog 1.0", usage = "\n\n    %prog --settings bia_pipeline.config [--target_task TASK] [more_options]")
    
                                
    #
    #   general options: verbosity / logging
    #
    parser.add_option("-v", "--verbose", dest = "verbose",
                      action="count", 
                      help="Print more verbose messages for each additional verbose level.")
    parser.add_option("-L", "--log_file", dest="log_file",
                      metavar="FILE",
                      type="string",
                      help="Name and path of log file")


    #
    #   pipeline
    #
    parser.add_option("-s", "--settings", dest="pipeline_settings",
                        metavar="FILE",
                        type="string",
                        help="File containing all the settings for the analysis.")                  
    parser.add_option("-t", "--target_tasks", dest="target_tasks",
                        action="append",
                        metavar="JOBNAME",
                        type="string",
                        help="Target task(s) of pipeline.")
    parser.add_option("-j", "--jobs", dest="jobs",
                        metavar="N",
                        type="int",
                        help="Allow N jobs (commands) to run simultaneously.")
    parser.add_option("-n", "--just_print", dest="just_print",
                        action="store_true", 
                        help="Don't actually run any commands; just print the pipeline.")
    parser.add_option("--flowchart", dest="flowchart",
                        metavar="FILE",
                        type="string",
                        help="Don't actually run any commands; just print the pipeline "
                             "as a flowchart.")

    #
    #   Less common pipeline options
    #
    parser.add_option("--key_legend_in_graph", dest="key_legend_in_graph",
                        action="store_true",
                        help="Print out legend and key for dependency graph.")
    parser.add_option("--forced_tasks", dest="forced_tasks",
                        action="append",
                        metavar="JOBNAME",
                        type="string",
                        help="Pipeline task(s) which will be included even if they are up to date.")
    parser.add_option("--rebuild_mode", dest="rebuild_mode",
                        action="store_false", 
                        help="gnu_make_maximal_rebuild_mode")
    parser.add_option("--run_on_bcl_tile", dest="run_on_bcl_tile",
                        type="string",                        
                        help="Use only this tile when doing bcl2fastq conversion. For testing purposes.")
    
    parser.set_defaults(pipeline_settings=None, 
                        jobs=1, verbose=0, 
                        target_tasks=list(), forced_tasks=list(), 
                        just_print=False, key_legend_in_graph=False,
                        rebuild_mode=True, run_on_bcl_tile=None)
    

    # get help string
    f =StringIO.StringIO()
    parser.print_help(f)
    helpstr = f.getvalue()
    (options, remaining_args) = parser.parse_args()


    #vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
    #                                             #
    #   Change this if necessary                  #
    #                                             #
    #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    #
    #   Add names of mandatory options,
    #       strings corresponding to the "dest" parameter
    #       in the options defined above
    #
    mandatory_options = ['pipeline_settings']

    #vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
    #                                             #
    #   Change this if necessary                  #
    #                                             #
    #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


    def check_mandatory_options (options, mandatory_options, helpstr):
        """
        Check if specified mandatory options have been defined
        """
        missing_options = []
        for o in mandatory_options:
            if not getattr(options, o):
                missing_options.append("--" + o)
    
        if not len(missing_options):
            return
    
        raise Exception("Missing mandatory parameter%s: %s.\n\n%s\n\n" %
                        ("s" if len(missing_options) > 1 else "",
                         ", ".join(missing_options),
                         helpstr))
        
    check_mandatory_options(options, mandatory_options, helpstr)
            
    



#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   Logger


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888



if __name__ == '__main__':
    import logging
    import logging.handlers

    MESSAGE = 15
    logging.addLevelName(MESSAGE, "MESSAGE")

    def setup_std_logging (logger, log_file, verbose):
        """
        set up logging using programme options
        """
        class debug_filter(logging.Filter):
            """
            Ignore INFO messages
            """
            def filter(self, record):
                return logging.INFO != record.levelno

        class NullHandler(logging.Handler):
            """
            for when there is no logging
            """
            def emit(self, record):
                pass

        # We are interesting in all messages
        logger.setLevel(logging.DEBUG)
        has_handler = False

        # log to file if that is specified
        if log_file:
            handler = logging.FileHandler(log_file, delay=False)
            handler.setFormatter(logging.Formatter("%(asctime)s - %(name)s - %(levelname)6s - %(message)s"))
            handler.setLevel(MESSAGE)
            logger.addHandler(handler)
            has_handler = True

        # log to stderr if verbose
        if verbose:
            stderrhandler = logging.StreamHandler(sys.stderr)
            stderrhandler.setFormatter(logging.Formatter("    %(message)s"))
            stderrhandler.setLevel(logging.DEBUG)
            if log_file:
                stderrhandler.addFilter(debug_filter())
            logger.addHandler(stderrhandler)
            has_handler = True

        # no logging
        if not has_handler:
            logger.addHandler(NullHandler())


    #
    #   set up log
    #
    module_name = "bia_pipeline"
    logger = logging.getLogger(module_name)
    setup_std_logging(logger, options.log_file, options.verbose)

    #
    #   Allow logging across Ruffus pipeline
    #
    def get_logger (logger_name, args):
        return logger

    from ruffus.proxy_logger import *
    (logger_proxy,
     logging_mutex) = make_shared_logger_and_proxy (get_logger,
                                                    module_name,
                                                    {})



    


    #vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
    #                                             #
    #   Get pipeline settings from a config file  #
    #                                             #
    #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

	if not os.path.exists(options.pipeline_settings): 
		raise Exception('Provided settings file [%s] does not exist or cannot be read.' % options.pipeline_settings)

    import ConfigParser
    config = ConfigParser.ConfigParser()
    config.read(options.pipeline_settings)



    # Should dockerized execution be used?
    dockerize = True
    try:
        docker_bin = config.get('Docker','docker-binary')
        logger.info('Found docker-binary setting. Using dockerized execution mode.')
    except ConfigParser.NoOptionError:
        logger.info('Docker-binary setting is missing. Using regular execution mode.')
        dockerize = False
 

    # Get the pipeline input
    run_folder = None
    input_fastqs = None
    input_bams = None
    try:
		run_folder = config.get('Inputs','run-directory') 
		logger.info('Found run-directory setting. Starting from bcl2fastq conversion.')
        # check presence of the run folder, and sample sheet file
        if not os.path.exists(run_folder) or not os.path.exists(os.path.join(run_folder,'SampleSheet.csv')):
            raise Exception("Missing sample sheet file: %s.\n" % os.path.join(run_folder,'SampleSheet.csv'))
    except ConfigParser.NoOptionError:
        try:
            input_fastqs = os.path.join(runs_scratch_dir if dockerize else '', config.get('Inputs','input-fastqs'))
            input_fastqs_resolved = glob.glob(input_fastqs)
            if len(input_fastqs_resolved) < 2:
                raise Exception("Missing input FASTQs. Found %s FASTQ files in [%s].\n" % (len(input_fastqs_resolved), config.get('Inputs','input-fastqs')))
            logger.info('Found %s input FASTQs. Starting from read trimming.' % len(input_fastqs_resolved))
        except ConfigParser.NoOptionError:
            try:
                input_bams = os.path.join(runs_scratch_dir if dockerize else '', config.get('Inputs','input-bams'))
                input_bams_resolved = glob.glob(input_bams)
                if len(input_bams_resolved) < 1:
                    raise Exception("Missing input BAMs. Found %s BAM files in [%s].\n" % (len(input_bams_resolved), config.get('Inputs','input-bams')))
                logger.info('Found %s input BAMs. Starting from indexing BAMs.' % len(input_bams_resolved))
            except ConfigParser.NoOptionError:
	        raise Exception('Found no valid input setting in [%s]. Please provide one of [run_directory|input-fastqs|input-bams] in the pipeline settings file.' % options.pipeline_settings)


 


    # Root paths
    
    reference_root = config.get('Paths','reference-root')
    scratch_root = os.getcwd()
    try:
        scratch_root   = config.get('Paths','scratch-root')
    except ConfigParser.NoOptionError:
        logger.info('Scratch-root setting is missing. Using current directory: %s' % scratch_root)
    
    run_id = os.path.basename(run_folder) if run_folder != None else os.path.basename(scratch_root)
    runs_scratch_dir = os.path.join(scratch-root, run_id) if run_folder != None else scratch_root
    logger.info('Run\'s scratch directory: %s' % runs_scratch_dir)
      
    # optional results and fastq archive dirs  
    results_archive = None
    try:
        results_archive = config.get('Paths','results-archive')
    except ConfigParser.NoOptionError:
        logger.info('No results-archive provided. Results will not be archived outside of the run\'s scratch directory.')
    
    fastq_archive = None
    try:
        fastq_archive = config.get('Paths','fastq-archive')
    except ConfigParser.NoOptionError:
        logger.info('No fastq-archive provided. Fastq files will not be archived outside of the run\'s scratch directory.')

    
    # optional /tmp dir
    tmp_dir = None
    try:
        tmp_dir = config.get('Paths','tmp-dir')
    except ConfigParser.NoOptionError:
        logger.info('No tmp-dir provided. %s\'s /tmp will be used.' % ('Container' if dockerize else 'Execution host'))
    

    if dockerize:
        # Docker args
        docker_args = config.get('Docker', 'docker-args')
        docker_args += " -v " + ":".join([run_folder, run_folder,"ro"])
        docker_args += " -v " + ":".join([reference_root,reference_root,"ro"])
    
        # Mount archive dirs as files from them are read (linked fastqs, gvcfs). 
        # Archiving is not performed by docker, so no write access should be needed.
        if fastq_archive != None:
            docker_args += " -v " + ":".join([fastq_archive,fastq_archive,"ro"])
        if results_archive != None:
            docker_args += " -v " + ":".join([results_archive,results_archive,"ro"])
    
        # Tmp, if should be different than the default  
        if tmp_dir != None: 
            docker_args += " -v " + ":".join([tmp_dir,tmp_dir,"rw"])
            
        docker_args += " -v " + ":".join([runs_scratch_dir,runs_scratch_dir,"rw"])
        docker_args += " -w " + runs_scratch_dir
        docker = " ".join([docker_bin, docker_args]) 
    
    # set the default value if the tmp-dir was unset
    tmp_dir = "/tmp" if tmp_dir==None else tmp_dir
     
    
    # reference files
    reference = os.path.join(reference_root, config.get('Resources','reference-genome'))    
    adapters = os.path.join(reference_root, config.get('Resources', 'adapters-fasta'))
    
    # tools
    bcl2fastq = config.get('Tools','bcl2fastq')
    trimmomatic = config.get('Tools', 'trimmomatic') 
	spades = config.get('Tools', 'spades') 

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   Logger


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888



if __name__ == '__main__':
    import logging
    import logging.handlers

    MESSAGE = 15
    logging.addLevelName(MESSAGE, "MESSAGE")

    def setup_std_logging (logger, log_file, verbose):
        """
        set up logging using programme options
        """
        class debug_filter(logging.Filter):
            """
            Ignore INFO messages
            """
            def filter(self, record):
                return logging.INFO != record.levelno

        class NullHandler(logging.Handler):
            """
            for when there is no logging
            """
            def emit(self, record):
                pass

        # We are interesting in all messages
        logger.setLevel(logging.DEBUG)
        has_handler = False

        # log to file if that is specified
        if log_file:
            handler = logging.FileHandler(log_file, delay=False)
            handler.setFormatter(logging.Formatter("%(asctime)s - %(name)s - %(levelname)6s - %(message)s"))
            handler.setLevel(MESSAGE)
            logger.addHandler(handler)
            has_handler = True

        # log to stderr if verbose
        if verbose:
            stderrhandler = logging.StreamHandler(sys.stderr)
            stderrhandler.setFormatter(logging.Formatter("    %(message)s"))
            stderrhandler.setLevel(logging.DEBUG)
            if log_file:
                stderrhandler.addFilter(debug_filter())
            logger.addHandler(stderrhandler)
            has_handler = True

        # no logging
        if not has_handler:
            logger.addHandler(NullHandler())


    #
    #   set up log
    #
    module_name = "exome"
    logger = logging.getLogger(module_name)
    setup_std_logging(logger, options.log_file, options.verbose)

    #
    #   Allow logging across Ruffus pipeline
    #
    def get_logger (logger_name, args):
        return logger

    from ruffus.proxy_logger import *
    (logger_proxy,
     logging_mutex) = make_shared_logger_and_proxy (get_logger,
                                                    module_name,
                                                    {})






#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#  Common functions 


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

import drmaa
drmaa_session = drmaa.Session()
drmaa_session.initialize()

from ruffus.drmaa_wrapper import run_job, error_drmaa_job

   
"""
cmd is given in a form:
    
    command {args}
    interpreter {interpreter_args} command {atgs}

The strings {args} and {interpreter_args} are replaced with args and interpreter_args values.
Examples of correct commands:
    cmd = "bcl2fastq {args}"
    cmd = "bcl2fastq -p param {args}"
    cmd = "bcl2fastq -p2 param {args} -p2 param2"
    cmd = "java {interpreter_args} -jar myjarfile.jar {args} -extras extra_arg
    cmd = "java -XmX4G {interpreter_args} -jar myjarfile.jar {args} -extras extra_arg
"""
def run_cmd(cmd, args, dockerize, interpreter_args=None, run_locally=True,
            cpus=1, mem_per_cpu=1024, walltime='24:00:00', 
            retain_job_scripts = True, job_script_dir = os.path.join(runs_scratch_dir, "drmaa")):
    
    full_cmd = ("{docker} "+cmd).format(docker = docker if dockerize else "",
                                        args=args, 
                                        interpreter_args = interpreter_args if interpreter_args!=None else "")

    stdout, stderr = "", ""
    job_options = "--account={account} \
		   --ntasks=1 \
                   --cpus-per-task={cpus} \
                   --mem-per-cpu={mem} \
                   --time={time} \
                  ".format(account=tsd_account, cpus=cpus, mem=int(1.2*mem_per_cpu), time=walltime)
                   
    try:
        stdout, stderr = run_job(full_cmd.strip(), 
                                 job_other_options=job_options,
                                 run_locally = run_locally, 
                                 retain_job_scripts = retain_job_scripts, job_script_directory = job_script_dir,
                                 logger=logger, working_directory=os.getcwd(),
                                 drmaa_session = drmaa_session)
    except error_drmaa_job as err:
        raise Exception("\n".join(map(str, ["Failed to run:", cmd, err, stdout, stderr])))


def get_sample_ids():
    """ Provides meaningful result only after HaplotypeCaller step"""
    files = glob.glob(os.path.join(runs_scratch_dir,'*','*.gvcf'))
    return [ os.path.splitext(os.path.basename(f))[0] for f in files ]

def get_num_files():
    """ Provides meaningful result only after HaplotypeCaller step"""
    return len(get_sample_ids())
    


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   Pipeline


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

from ruffus import *


#
#
# Prepare FASTQ
# 

@active_if(run_folder != None)
@follows(mkdir(runs_scratch_dir), mkdir(os.path.join(runs_scratch_dir,'fastqs')))
@files(run_folder, os.path.join(runs_scratch_dir,'fastqs','completed'))
@posttask(touch_file(os.path.join(runs_scratch_dir,'fastqs','completed')))
def bcl2fastq_conversion(run_directory, completed_flag):
    """ Run bcl2fastq conversion and create fastq files in the run directory"""
    out_dir = os.path.join(runs_scratch_dir,'fastqs')
    interop_dir = os.path.join(out_dir,'InterOp')

    # r, w, d, and p specify numbers of threads to be used for each of the concurrent subtasks of the conversion (see bcl2fastq manual) 
    args = "-R {indir} -o {outdir} --interop-dir={interopdir} -r1 -w1 -d2 -p4 \
            ".format(indir=run_directory, outdir=out_dir, interopdir=interop_dir)
    if options.run_on_bcl_tile != None:
        args += " --tiles %s" % options.run_on_bcl_tile
        
    run_cmd(bcl2fastq, args, dockerize=dockerize, cpus=8, mem_per_cpu=2048)
    


@active_if(run_folder != None and fastq_archive != None)
@transform(bcl2fastq_conversion, formatter(".+/(?P<RUN_ID>[^/]+)/fastqs/completed"), str(fastq_archive)+"/{RUN_ID[0]}")
def archive_fastqs(completed_flag, archive_dir):
    """ Archive fastqs """    
    fq_dir = os.path.dirname(completed_flag)

# uncomment if archive should not be overwritten (risk of creating many archives of the same run)
#    if os.path.exists(archive_dir):
#	import time
#	archive_dir += "_archived_"+str(time.strftime("%Y%m%d_%H%M%S"))

    import shutil
    shutil.move(fq_dir, archive_dir)
    os.mkdir(fq_dir)
    for f in glob.glob(os.path.join(archive_dir,"*.fastq.gz")):
        os.symlink(f, os.path.join(fq_dir,os.path.basename(f)))


#
# Prepare directory for every sample and link the input fastq files
# Expected format:
#    /path/to/file/[SAMPLE_ID]_S[1-9]\d?_L\d\d\d_R[12]_001.fastq.gz
# SAMPLE_ID can contain all signs except path delimiter, i.e. "\"
#
@active_if(run_folder != None or input_fastqs != None)
@jobs_limit(1)    # to avoid problems with simultanous creation of the same sample dir
@follows(archive_fastqs)
@subdivide(os.path.join(runs_scratch_dir,'fastqs','*.fastq.gz') if run_folder != None else input_fastqs,
           formatter('(?P<PATH>.+)/(?P<SAMPLE_ID>[^/]+)_S[1-9]\d?_L\d\d\d_R[12]_001\.fastq\.gz$'), 
           '{subpath[0][1]}/{SAMPLE_ID[0]}/{basename[0]}{ext[0]}')
def link_fastqs(fastq_in, fastq_out):
    """Make working directory for every sample and link fastq files in"""
    if not os.path.exists(os.path.dirname(fastq_out)):
        os.mkdir(os.path.dirname(fastq_out))
    if not os.path.exists(fastq_out):
        os.symlink(fastq_in, fastq_out) 

    
    
#
# Input FASTQ filenames are expected to have following format:
#    [SAMPLE_ID]_[S_NUM]_[LANE_ID]_[R1|R2]_001.fastq.gz
# In this step, the two FASTQ files matching on the [SAMPLE_ID]_[S_ID]_[LANE_ID] will be trimmed together (R1 and R2). 
# The output will be written to two FASTQ files
#    [SAMPLE_ID]_[LANE_ID].fq1.gz
#    [SAMPLE_ID]_[LANE_ID].fq2.gz
# SAMPLE_ID can contain all signs except path delimiter, i.e. "\"
#
@active_if(run_folder != None or input_fastqs != None)
@collate(link_fastqs, regex(r'(.+)/([^/]+)_S[1-9]\d?_(L\d\d\d)_R[12]_001\.fastq\.gz$'),  r'\1/\2_\3.fq1.gz')
def trim_reads(inputs, output):
    outfq1 = output
    outfq2 = output.replace('fq1.gz','fq2.gz')
    unpaired = [outfq1.replace('fq1.gz','fq1_unpaired.gz'), outfq2.replace('fq2.gz','fq2_unpaired.gz')]               
    # logfile = output.replace('fq1.gz','trimmomatic.log')
    # -trimlog {log} \
    # log=logfile
    args = "PE -phred33 -threads 1 \
            {in1} {in2} {out1} {unpaired1} {out2} {unpaired2} \
            MINLEN:36 \
            ILLUMINACLIP:{adapter}:2:30:10 \
            LEADING:3 \
            TRAILING:3 \
            SLIDINGWINDOW:4:15".format(in1=inputs[0], in2=inputs[1],
                                       out1=outfq1, out2=outfq2,
                                       unpaired1=unpaired[0], unpaired2=unpaired[1],
                                       adapter=adapters)
    max_mem = 2048
    run_cmd(trimmomatic, args, interpreter_args="-Xmx"+str(max_mem)+"m", 
            dockerize=dockerize, cpus=1, mem_per_cpu=max_mem)


#
#
# Assemble the reads 
# 


def clean_trimmed_fastqs():
    """ Remove the trimmed fastq files. Links to original fastqs are kept """
    for f in glob.glob(os.path.join(runs_scratch_dir,'*','*.fq[12]*.gz')):
        os.remove(f)

#
# FASTQ filenames are expected to have following format:
#    [SAMPLE_ID]_[LANE_ID].fq[1|2].gz
# In this step, the fq1 file coming from trim_reads is matched with the fq2 file and assembled together. 
# The output will be written to SAMPLE_ID directory:
#    [SAMPLE_ID]/
#
@jobs_limit(8)
#@posttask(clean_trimmed_fastqs)
@collate(trim_reads, regex(r'(.+)/([^/]+)_S[1-9]\d?_L\d\d\d_R[12]_001\.fastq\.gz$'),  r'\2/contigs.fasta')
def assemble_reads(fastqs, contigs):
    threads = 4
    
    out_dir=os.path.dirname(contigs)
    
    args = "-1 {fq1} -2 {fq2} -o {out_dir} \
		    -m {mem} -t {threads} --carefull \
			".format(fq1=fastqs[0], fq2=fastqs[1], out_dir=out_dir, threads=threads)

    run_cmd(spades, args, cpus=threads, mem_per_cpu=8192/threads)
    
    
    
    

#
#
# QC the FASTQ files
#

@transform(link_fastqs, suffix(".gz"), '.html')
def qc_raw_reads(input_fastq, report):
    """ Generate FastQC report for raw FASTQs """
    produce_fastqc_report(input_fastq, report)

@transform(trim_reads, suffix(".gz"), '.html')
def qc_trimmed_reads(input_fastq, report):
    """ Generate FastQC report for raw FASTQs """
    produce_fastqc_report(input_fastq, report)

@follows(qc_raw_reads, qc_trimmed_reads)
def qc_reads():
	pass

#
#
# QC the assemblies
#

@merge(assemble_reads, os.path.join(runs_scratch_dir, 'qc', run_id))
def qc_assemblies(contigs, report_dir):
	args = "-o {%s}" % report_dir + " ".join(contigs)
    run_cmd(quast, args)




#def archive_results():
    ## if optional results_archive was not provided - do nothing
    #if results_archive == None: return
    #arch_path = os.path.join(results_archive, run_id)
    #if not os.path.exists(arch_path): 
        #os.mkdir(arch_path)
        
    #run_cmd("cp %s/*/*.gatk.bam %s" % (runs_scratch_dir,arch_path), "", run_locally=True)
    #run_cmd("cp %s/*/*.gatk.bam.gene_coverage* %s" % (runs_scratch_dir,arch_path), "", run_locally=True)
    #run_cmd("cp %s/*/*.exome.vcf %s" % (runs_scratch_dir,arch_path), "", run_locally=True)
    #run_cmd("cp %s/*.multisample.gvcf %s" % (runs_scratch_dir, results_archive),
            #"", run_locally=True)
    #run_cmd("cp -r %s/qc %s" % (runs_scratch_dir,arch_path), "", run_locally=True)


#def cleanup_files():
    #run_cmd("rm -rf {dir}/*/*.recal_data.csv {dir}/*/*.realign* {dir}/*/*.dedup* \
            #{dir}/*.multisample.indel.model* {dir}/*.multisample.snp.model* \
            #{dir}/*/*.log {dir}/*.multisample.recalibratedSNPS.rawIndels.vcf* \
            #{dir}/*.multisample.recalibrated.vcf* \
            #".format(dir=runs_scratch_dir), "", run_locally=True)


#@posttask(archive_results, cleanup_files)
@follows(qc_reads, qc_assemblies)
def complete_run():
    pass





#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   Main logic


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
if __name__ == '__main__':
    if options.just_print:
        pipeline_printout(sys.stdout, options.target_tasks, options.forced_tasks,
                            gnu_make_maximal_rebuild_mode = options.rebuild_mode,
                            verbose=options.verbose, #verbose_abbreviated_path=0,
                            checksum_level = 0)

    elif options.flowchart:
        pipeline_printout_graph (   open(options.flowchart, "w"),
                                    # use flowchart file name extension to decide flowchart format
                                    #   e.g. svg, jpg etc.
                                    os.path.splitext(options.flowchart)[1][1:],
                                    options.target_tasks,
                                    options.forced_tasks,
                                        gnu_make_maximal_rebuild_mode = options.rebuild_mode,
                                    no_key_legend   = not options.key_legend_in_graph)
    else:        
        pipeline_run(options.target_tasks, options.forced_tasks,
                            multithread     = options.jobs,
                            logger          = stderr_logger,
                            verbose         = options.verbose,
                            gnu_make_maximal_rebuild_mode = options.rebuild_mode,
                            checksum_level  = 0)
    
        
    drmaa_session.exit()
    
