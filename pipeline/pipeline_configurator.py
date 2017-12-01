#
#
#
# Pipeline configurator
#

import os, json

class PipelineConfig:


    # Attributes:
    #
    #  - reference_root
    #  - scratch_root
    #  - results_archive
    #  - fastq_archive
    #  - tmp_dir
    


    # Tool paths
    # - docker_bin
    

    def __init__(self):
        
        self.drmaa_session = None
        
        self.dockerize = True
        
        self.logger = None
        
        self.bcl2fastq   = None
        self.bbmerge     = None
        self.trimmomatic = None
        self.fastqc      = None
        self.spades      = None
        self.quast       = None
        self.bwa         = None
        self.samtools    = None
        self.freebayes   = None

        # run settings
        self.reference_root = None
        self.scratch_root = None
        self.run_folder   = None
        self.input_fastqs = None
        self.run_id       = None
        self.runs_scratch_dir = None
        self.tmp_dir      = None

        self.adapters  = None
        self.reference = None

        self.target_tasks    = []
        self.log_file        = None
        self.num_jobs        = 1
        self.verbosity_level = 0
        self.dry_run         = False
        self.rebuild_mode    = False
        self.run_on_bcl_tile = None


    def set_logger(self, logger):
        self.logger = logger
   
    def set_runfolder(self, runfolder):
        if runfolder != None and \
            os.path.exists(runfolder) and \
            os.path.exists(os.path.join(runfolder,'SampleSheet.csv')):

            self.run_folder = runfolder
            
        else:
            raise Exception("Incorrect runfolder\'s path [%s] or missing SampleSheet file." % run_folder)
   

    def set_input_fastqs(self, fastqs):
        self.input_fastqs = fastqs
    
    
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
        try:
            self.docker_bin = config.get('Docker','docker-binary')
            self.logger.info('Found docker-binary setting. Using dockerized execution mode.')
        except ConfigParser.NoOptionError:
            self.logger.info('Docker-binary setting is missing. Using regular execution mode.')
            self.dockerize = False
    
        # Get the pipeline input    
        if self.run_folder is None:
            try:
                self.run_folder = config.get('Inputs','run-directory') 
                logger.info('Found run-folder setting. Starting from bcl2fastq conversion of %s.' % self.run_folder)
                # check presence of the run folder, and sample sheet file
                if not os.path.exists(self.run_folder) or not os.path.exists(os.path.join(self.run_folder,'SampleSheet.csv')):
                    raise Exception("Missing sample sheet file: %s.\n" % os.path.join(self.run_folder,'SampleSheet.csv'))
            except ConfigParser.NoOptionError:
                try:
                    input_fastqs = os.path.join(self.runs_scratch_dir if self.dockerize else '', config.get('Inputs','input-fastqs'))
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
        
        self.run_id = os.path.basename(self.run_folder) if self.run_folder != None else os.path.basename(self.scratch_root)
        
        #
        # TODO
        # needs to be updated on update of settings
        #
        self.runs_scratch_dir = os.path.join(self.scratch_root, self.run_id) if self.run_folder != None else self.scratch_root
        self.logger.info('Run\'s scratch directory: %s' % self.runs_scratch_dir)
          
        # optional results and fastq archive dirs  
        self.results_archive = None
        try:
            self.results_archive = config.get('Paths','results-archive')
        except ConfigParser.NoOptionError:
            self.logger.info('No results-archive provided. Results will not be archived outside of the run\'s scratch directory.')
        
        self.fastq_archive = None
        try:
            self.fastq_archive = config.get('Paths','fastq-archive')
        except ConfigParser.NoOptionError:
            self.logger.info('No fastq-archive provided. Fastq files will not be archived outside of the run\'s scratch directory.')
    
        
        # optional /tmp dir
        try:
            self.tmp_dir = config.get('Paths','tmp-dir')
        except ConfigParser.NoOptionError:
            self.logger.info('No tmp-dir provided. %s\'s /tmp will be used.' % ('Container' if self.dockerize else 'Execution host'))
        
    
        #if self.dockerize:
            ## Docker args
            #docker_args = config.get('Docker', 'docker-args')
            #docker_args += " -v " + ":".join([run_folder, run_folder,"ro"])
            #docker_args += " -v " + ":".join([self.reference_root,self.reference_root,"ro"])
        
            ## Mount archive dirs as files from them are read (linked fastqs, gvcfs). 
            ## Archiving is not performed by docker, so no write access should be needed.
            #if fastq_archive != None:
                #self.docker_args += " -v " + ":".join([fastq_archive,fastq_archive,"ro"])
            #if results_archive != None:
                #self.docker_args += " -v " + ":".join([results_archive,results_archive,"ro"])
        
            ## Tmp, if should be different than the default  
            #if self.tmp_dir != None: 
                #self.docker_args += " -v " + ":".join([self.tmp_dir,self.tmp_dir,"rw"])
                
            #self.docker_args += " -v " + ":".join([self.runs_scratch_dir,self.runs_scratch_dir,"rw"])
            #self.docker_args += " -w " + self.runs_scratch_dir
            #self.docker = " ".join([self.docker_bin, self.docker_args]) 
        
        # set the default value if the tmp-dir was unset
        self.tmp_dir = "/tmp" if self.tmp_dir==None else self.tmp_dir
         
        
        # reference files
        self.reference = os.path.join(self.reference_root, config.get('Resources','reference-genome'))    
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


        

    def load_settings_from_JSON(self, settings):
        """ Update setting with JSON dictionary """
        
        d = json.loads(settings, parse_int=int)
        if d.has_key('target_tasks'):
            self.target_tasks = d['target_tasks']
        if d.has_key('num_jobs'):
            self.num_jobs = int(d['num_jobs'])
    

    def is_runnable(self):
        """ returns true if all required settings are set """
        
        if len(target_tasks) < 1: 
            return False
        # check task names?
        
        if self.run_folder is None or \
            not os.path.exists(self.run_folder) or \
            not os.path.exists(os.path.join(self.run_folder, self.run_id, 'SampleSheet.csv')):
            return False
            
        return True


    def write_config(self, config_file):
        """ Save the config in a file """
        
        # write root paths
        
        # write reference data
        
        # write tool paths
        
        pass
        
        
