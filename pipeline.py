#
#
#
# 1) Contacts the database to learn which tasks should be run on which samples
#   
#   a dictionary of (sample_ID: settings) is returned
#   samples are grouped by setttings values (eg. mr_assembly_qc)
#  
# 2) for each group of samples / pipeline run a PipelineConfig object is created
#
# 3) sequenctailly pipeline is executed with the config objects
#  
#
#
#






# 8888888888888888888888888888888888888888
#
#           C L I   a r g s
#
# 8888888888888888888888888888888888888888



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
    
def parse_cli_args():
    
    from optparse import OptionParser
    import StringIO

    parser = OptionParser(version="%prog 1.0", usage = "\n\n    %prog --settings bia_pipeline.config --run_folder RUN_FOLDER [--target_task TASK] [more_options]")
    
                                
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
    #   pipeline config and control
    #
    parser.add_option("-s", "--settings", dest="pipeline_settings",
                        metavar="FILE",
                        type="string",
                        help="File containing all the settings for the analysis.")                  
    parser.add_option("-r", "--run_folder", dest="run_folder",
                        metavar="FILE",
                        type="string",
                        help="Directory with run's raw data as produced by Illumina sequencer")                  
    parser.add_option("-t", "--target_tasks", dest="target_tasks",
                        action="append",
                        metavar="JOBNAME",
                        type="string",
                        help="Target task(s) of pipeline.")
#    parser.add_option("-j", "--jobs", dest="jobs",
#                        metavar="N",
#                        type="int",
#                        help="Allow N jobs (commands) to run simultaneously.")
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
    f = StringIO.StringIO()
    parser.print_help(f)
    helpstr = f.getvalue()
    (options, remaining_args) = parser.parse_args()


    mandatory_options = ['pipeline_settings', 'run_folder']
    check_mandatory_options(options, mandatory_options, helpstr)
 
    return options






# 8888888888888888888888888888888888888888
#
#         L o g g i n g
#
# 8888888888888888888888888888888888888888

import sys
import logging
import logging.handlers


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


def setup_logging(log_file, verbosity):

    MESSAGE = 15
    logging.addLevelName(MESSAGE, "MESSAGE")
    
    module_name = "bia_pipeline"
    logger = logging.getLogger(module_name)
    setup_std_logging(logger, log_file, verbosity)

    return logger

    
   


# 8888888888888888888888888888888888888888
#
#   Per-sample config of the pipeline
#
# 8888888888888888888888888888888888888888



class SampleTable:
    
    def get_sample_configs(self, sample_ids):
        """ returns a dict with sample_ID : string. 
        The string is a JSON key:value list of run settings for the sample.
        Among supported keys are: target_task (string), num_jobs (integer), ..."""
        
        return {s:'{"target_tasks":["assemble_merged"]}' for s in sample_ids}



def read_sample_ids_from_samplesheet():
    
    sample_ids = []
    with open(os.path.join(options.run_folder, 'SampleSheet.csv')) as ss:
        
        lines = [l.strip() for l in ss.readlines()]
        header_idx = lines.index('[Data]')+1
        
        # find index of sample ID (most likely 0)
        header = lines[header_idx].split(',')
        sample_id_index = header.index('Sample_ID')
                
        for l in lines[header_idx+1:len(lines)]:
            lsplit = l.split(',')
            if len(lsplit) == len(header):
                sample_ids.append(lsplit[sample_id_index])
        
    return sample_ids    
    
    
   
# 88888888888888888888888888888888888888888
#
#              M A I N
#
# 88888888888888888888888888888888888888888


if __name__ == '__main__':



    options = parse_cli_args()


    logger = setup_logging(options.log_file, options.verbose)
    
    #
    #   Allow logging across Ruffus pipeline
    #
    def get_logger(logger_name, args):
        return logger

    from ruffus.proxy_logger import *
    (logger_proxy, logging_mutex) = \
        make_shared_logger_and_proxy(get_logger, logger.name, {})
        


    sample_ids = read_sample_ids_from_samplesheet()

    st = SampleTable()
    cfgs = st.get_sample_configs(sample_ids)
    
    # group by configs
    cfgs_groups = {}
    for k, v in cfgs.iteritems():
        cfgs_groups.setdefault(v, []).append(k)

    
    import pipeline.global_vars
    
    import drmaa
    pipeline.global_vars.drmaa_session = drmaa.Session()
    pipeline.global_vars.drmaa_session.initialize()

    from ruffus import pipeline_printout, pipeline_printout_graph, pipeline_run
    from pipeline.pipeline_configurator import PipelineConfig
    
    
    for (cfg_group_idx, cfg_group) in enumerate(cfgs_groups.keys()):
       
        pipeline.global_vars.cfg = PipelineConfig() 
        pipeline.global_vars.cfg.set_logger(logger)
        pipeline.global_vars.cfg.set_runfolder(options.run_folder)
        pipeline.global_vars.cfg.load_settings_from_file(options.pipeline_settings)
        pipeline.global_vars.cfg.load_settings_from_JSON(cfg_group) # cfg_group is the group identifier and JSON string at the same time
    

        #
        #
        # task definitions should use e.g. cfg.reference
        #
        #

        from pipeline.tasks import *


        if options.just_print:
            
            print "SAMPLE_GROUP_" + str(cfg_group_idx) + ":"
            
            pipeline_printout(sys.stdout, cfg.target_tasks, options.forced_tasks,
                            gnu_make_maximal_rebuild_mode = options.rebuild_mode,
                            verbose=options.verbose, verbose_abbreviated_path=0,
                            checksum_level = 0)

        elif options.flowchart:
                        
            pipeline_printout_graph ( open("SAMPLE_GROUP_" + \
                                            str(cfg_group_idx) + "_" + \
                                            options.flowchart, "w"),
                                    # use flowchart file name extension to decide flowchart format
                                    #   e.g. svg, jpg etc.
                                    os.path.splitext(options.flowchart)[1][1:],
                                    cfg.target_tasks, options.forced_tasks,
                                    gnu_make_maximal_rebuild_mode = options.rebuild_mode,
                                    no_key_legend = not options.key_legend_in_graph)
        else:        
            pipeline_run(cfg.target_tasks, options.forced_tasks,
                            multithread     = cfg.num_jobs,
                            logger          = logger,
                            verbose         = options.verbose,
                            gnu_make_maximal_rebuild_mode = options.rebuild_mode,
                            checksum_level  = 0)
    
       
    global_vars.drmaa_session.exit()

       
