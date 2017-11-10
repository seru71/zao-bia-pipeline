#
#
#
# Pipeline configurator
#

 
 
class PipelineConfig:


    # Attributes:
    #
    #  - reference_root
    #  - scratch_root
    #  - results_archive
    #  - fastq_archive
    #  - tmp_dir
    
    #  - adapters 
    #  - reference


    # Tool paths
    # - docker_bin
    
    bcl2fastq   = None
    bbmerge     = None
    trimmomatic = None
    fastqc      = None
    spades      = None
    quast       = None
    bwa         = None
    samtools    = None
    freebayes   = None

    # run settings
    run_folder   = None
    input_fastqs = None
    run_id       = None
    
    target_tasks    = []
    log_file        = None
    num_threads     = 1
    verbosity_level = 0
    dry_run         = False
    rebuild_mode    = False


    def __init__(self, logger):
        
        self.logger = logger
   
   
    def set_runfolder(self, runfolder):
        if runfolder != None and \
            os.path.exists(runfolder) and \
            os.path.exists(os.path.join(options.runfolder,'SampleSheet.csv')):

            self.run_folder = runfolder
            
        else:
            raise Exception("Incorrect runfolder\'s path [%s] or missing SampleSheet file." % run_folder)
    
    
    
    def load_settings_from_file(self, cfg_file):
        
        #
        #
        # TODO
        # Missing settings should not cause exceptions
        #
        #
        #
        
        
        
        """ preset settings using config file """

        if not os.path.exists(cfg_file): 
            raise Exception('Provided config file [%s] does not exist or cannot be read.' % cfg_file)

        import ConfigParser
        config = ConfigParser.ConfigParser()
        config.read(cfg_file)
        
        
        # Should dockerized execution be used?
        dockerize = True
        try:
            self.docker_bin = config.get('Docker','docker-binary')
            self.logger.info('Found docker-binary setting. Using dockerized execution mode.')
        except ConfigParser.NoOptionError:
            self.logger.info('Docker-binary setting is missing. Using regular execution mode.')
            dockerize = False
    
        # Get the pipeline input    
        if self.run_folder is None:
            try:
                self.run_folder = config.get('Inputs','run-directory') 
                logger.info('Found run-folder setting. Starting from bcl2fastq conversion of %s.' % self.run_folder)
                # check presence of the run folder, and sample sheet file
                if not os.path.exists(run_folder) or not os.path.exists(os.path.join(run_folder,'SampleSheet.csv')):
                    raise Exception("Missing sample sheet file: %s.\n" % os.path.join(run_folder,'SampleSheet.csv'))
            except ConfigParser.NoOptionError:
                try:
                    input_fastqs = os.path.join(runs_scratch_dir if dockerize else '', config.get('Inputs','input-fastqs'))
                    input_fastqs_resolved = glob.glob(input_fastqs)
                    if len(input_fastqs_resolved) < 2:
                        raise Exception("Missing input FASTQs. Found %s FASTQ files in [%s].\n" % (len(input_fastqs_resolved), config.get('Inputs','input-fastqs')))
                    self.logger.info('Found %s input FASTQs. Starting from read trimming.' % len(input_fastqs_resolved))
                except ConfigParser.NoOptionError:
                    
                    #
                    #
                    # TODO
                    # This shoudl not raise exeption as new setting can be added later
                    #
                    #
                    
                    
                    raise Exception('Found no valid input setting in [%s]. Please provide one of [run_directory|input-fastqs] in the pipeline settings file.' % options.pipeline_settings)
    
            
            
        
        
        
        
        
        
        
        
        
        
        self.reference_root = config.get('Paths','reference-root')
        self.scratch_root = os.getcwd()
        try:
            self.scratch_root   = config.get('Paths','scratch-root')
        except ConfigParser.NoOptionError:
            self.logger.info('Scratch-root setting is missing. Using current directory: %s' % self.scratch_root)
        
        self.run_id = os.path.basename(run_folder) if run_folder != None else os.path.basename(self.scratch_root)
        
        #
        # TODO
        # needs to be updated on update of settings
        #
        self.runs_scratch_dir = os.path.join(scratch_root, run_id) if self.run_folder != None else self.scratch_root
        self.logger.info('Run\'s scratch directory: %s' % self.runs_scratch_dir)
          
        # optional results and fastq archive dirs  
        self.results_archive = None
        try:
            rself.esults_archive = config.get('Paths','results-archive')
        except ConfigParser.NoOptionError:
            self.logger.info('No results-archive provided. Results will not be archived outside of the run\'s scratch directory.')
        
        self.fastq_archive = None
        try:
            self.fastq_archive = config.get('Paths','fastq-archive')
        except ConfigParser.NoOptionError:
            self.logger.info('No fastq-archive provided. Fastq files will not be archived outside of the run\'s scratch directory.')
    
        
        # optional /tmp dir
        self.tmp_dir = None
        try:
            tmp_dir = config.get('Paths','tmp-dir')
        except ConfigParser.NoOptionError:
            self.logger.info('No tmp-dir provided. %s\'s /tmp will be used.' % ('Container' if dockerize else 'Execution host'))
        
    
        #if dockerize:
            ## Docker args
            #docker_args = config.get('Docker', 'docker-args')
            #docker_args += " -v " + ":".join([run_folder, run_folder,"ro"])
            #docker_args += " -v " + ":".join([reference_root,reference_root,"ro"])
        
            ## Mount archive dirs as files from them are read (linked fastqs, gvcfs). 
            ## Archiving is not performed by docker, so no write access should be needed.
            #if fastq_archive != None:
                #docker_args += " -v " + ":".join([fastq_archive,fastq_archive,"ro"])
            #if results_archive != None:
                #docker_args += " -v " + ":".join([results_archive,results_archive,"ro"])
        
            ## Tmp, if should be different than the default  
            #if tmp_dir != None: 
                #docker_args += " -v " + ":".join([tmp_dir,tmp_dir,"rw"])
                
            #docker_args += " -v " + ":".join([runs_scratch_dir,runs_scratch_dir,"rw"])
            #docker_args += " -w " + runs_scratch_dir
            #docker = " ".join([docker_bin, docker_args]) 
        
        # set the default value if the tmp-dir was unset
        tmp_dir = "/tmp" if tmp_dir==None else tmp_dir
         
        
        # reference files
        self.reference = os.path.join(reference_root, config.get('Resources','reference-genome'))    
        self.adapters = config.get('Resources', 'adapters-fasta')
        
        # tools
        self.bcl2fastq = config.get('Tools','bcl2fastq')
        self.fastqc = config.get('Tools', 'fastqc')
        self.bbmerge = config.get('Tools', 'bbmerge') 
        self.trimmomatic = config.get('Tools', 'trimmomatic') 
        self.bwa = config.get('Tools', 'bwa') 
        self.samtools = config.get('Tools', 'samtools') 
        self.freebayes = config.get('Tools', 'freebayes') 
        self.spades = config.get('Tools', 'spades') 
        self.quast = config.get('Tools', 'quast')


        

    def load_setting_from_JSON(self, json_dictionary):
        """ Update setting with JSON dictionary """
        
        # update self


    def is_runnable(self):
        """ returns true if all required settings are set """
        
        if len(target_tasks) < 1: 
            return False
        # check task names?
        
        if runfolder is None or
            not os.path.exists(runfolder) or 
            not os.path.exists(os.path.join(run_folder, run_id, 'SampleSheet.csv'):
            return False
            
        return True


    def write_config(self, config_file):
        """ Save the config in a file """
        
        # write root paths
        
        # write reference data
        
        # write tool paths
        
        pass
        
        
        
