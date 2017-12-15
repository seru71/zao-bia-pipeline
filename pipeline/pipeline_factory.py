        
from ruffus import *
from pipeline.tasks import *
    
class ZaoPipelineFactory:    
    
    def __init__(self):
        pass
    
    
    def get_bacterial_assembly_pipeline(self, name, cfg):
    
        # assemble the pipeline 
        p = Pipeline(name=name)
        
        bcl2fastq_conversion_task = p.transform(bcl2fastq_conversion, 
                        cfg.run_folder, formatter(), 
                        os.path.join(cfg.runs_scratch_dir,'fastqs','completed'), 
                        cfg)\
                .follows(mkdir(cfg.runs_scratch_dir))\
                .follows(mkdir(os.path.join(cfg.runs_scratch_dir,'fastqs')))\
                .posttask(touch_file(os.path.join(cfg.runs_scratch_dir,'fastqs','completed')))\
                .posttask(lambda: log_task_progress(cfg, 'bcl2fastq_conversion', completed=True))


        archive_fastqs_task = p.transform(archive_fastqs, 
                                    bcl2fastq_conversion_task,
                                    formatter(".+/(?P<RUN_ID>[^/]+)/fastqs/completed"), 
                                    str(cfg.fastq_archive)+"/{RUN_ID[0]}/fastq")\
                               .active_if(cfg.fastq_archive != None)\
                               .posttask(lambda: log_task_progress(cfg, 'archive_fastqs', completed=True))
    
        
        link_fastqs_task = p.transform(link_fastqs,
                                    cfg.input_fastqs,
                                    formatter('(?P<PATH>.+)/(?P<SAMPLE_ID>[^/]+)_S[1-9]\d?_L\d\d\d_R[12]_001\.fastq\.gz$'),
                                    cfg.runs_scratch_dir+'/{SAMPLE_ID[0]}/{basename[0]}{ext[0]}')\
                            .follows(archive_fastqs_task)\
                            .jobs_limit(1)\
                            .posttask(lambda: log_task_progress(cfg, 'link_fastqs', completed=True))
                     
        
        
        
        # #######################################
        #
        #   R e a d   p r e p r o c e s s i n g                     
        #
        # ###################################
                                    
        #
        # just trim ...
        trim_reads_task = p.collate(trim_reads, 
                                    link_fastqs_task,
                                    regex(r'(.+)/([^/]+)_S[1-9]\d?_(L\d\d\d)_R[12]_001\.fastq\.gz$'),
                                    [r'\1/\2_\3_R1.fq.gz', r'\1/\2_\3_R2.fq.gz', 
                                    r'\1/\2_\3_R1_unpaired.fq.gz', r'\1/\2_\3_R2_unpaired.fq.gz'],
                                    cfg)\
                            .posttask(lambda: log_task_progress(cfg, 'trim_reads', completed=True))
   
        
        
        #
        # ... or merge & trim
        merge_reads_task = p.collate(merge_reads,
                                    link_fastqs_task,
                                    regex(r'(.+)/([^/]+)_S[1-9]\d?_(L\d\d\d)_R[12]_001\.fastq\.gz$'), 
                                    [r'\1/\2_merged.fq.gz', r'\1/\2_notmerged_R1.fq.gz', r'\1/\2_notmerged_R2.fq.gz'],
                                    cfg)\
                            .posttask(lambda: log_task_progress(cfg, 'merge_reads', completed=True))


        trim_notmerged_pairs_task = \
            p.transform(trim_notmerged_pairs, 
                        merge_reads_task,
                        formatter(None, '.+/(?P<PREFIX>[^/]+)\.fq\.gz$', '.+/(?P<PREFIX>[^/]+)\.fq\.gz$'), 
                        ['{path[1]}/{PREFIX[1]}.trimmed.fq.gz', '{path[2]}/{PREFIX[2]}.trimmed.fq.gz',
                        '{path[1]}/{PREFIX[1]}.unpaired.fq.gz', '{path[2]}/{PREFIX[2]}.unpaired.fq.gz'],
                        cfg)\
             .posttask(lambda: log_task_progress(cfg, 'trim_notmerged_pairs', completed=True))
                                                
        
        
        trim_merged_reads_task = \
            p.transform(trim_merged_reads, 
                        merge_reads_task,
                        suffix('_merged.fq.gz'), 
                        '_merged.trimmed.fq.gz',
                        cfg)\
            .posttask(lambda: log_task_progress(cfg, 'trim_merged_reads', completed=True))
   
   
        # #################
        #
        #   M a p p i n g
        #
        ###############
   
   
        map_trimmed_reads_task = \
            p.transform(map_trimmed_reads,
                        trim_reads_task,
                        formatter("(.+)/(?P<SAMPLE_ID>[^/]+)_L\d\d\d_R[12](_unpaired)?\.fq\.gz$"),
                        "{subpath[0][0]}/{subdir[0][0]}.bam", "{SAMPLE_ID[0]}",
			cfg)\
	     .posttask(lambda: log_task_progress(cfg, 'map_trimmed_reads', completed=True))
   
   
        # #################################
        #
        #   V a r i a n t   c a l l i n g
        #
        # #############################
        
      
        call_variants_on_trimmed_task = \
            p.transform(call_variants_on_trimmed,
                        map_trimmed_reads_task,
                        suffix(".bam"), 
                        ".fb.vcf",
			cfg)\
             .posttask(lambda: log_task_progress(cfg, 'call_variants_on_trimmed', completed=True))

   
        joincall_variants_on_triummed_task = \
            p.merge(jointcall_variants_on_trimmed,
                    map_trimmed_reads_task, 
                    os.path.join(cfg.runs_scratch_dir, "multisample.fb.vcf"),
                    cfg)\
             .posttask(lambda: log_task_progress(cfg, 'jointcall_variants_on_trimmed', completed=True))
   
   
        # ###################
        #
        #   A s s e m b l y
        #
        # ###############

        
        
#        assemble_trimmed_task = \
#            p.collate(assemble_trimmed,
#                      trim_reads_task,
#                      formatter(),
#                      '{subpath[0][0]}/{subdir[0][0]}_tra.fasta',
#		      cfg)\
#             .posttask(lambda: clean_assembly_dir('tra_assembly', cfg))\
#             .posttask(lambda: log_task_progress(cfg, 'assemble_trimmed', start=False))\
#             .jobs_limit(4)
                     
            
        assemble_merged_task = \
            p.collate(assemble_merged,
                      [trim_merged_reads, trim_notmerged_pairs], 
                      formatter(), 
                      '{subpath[0][0]}/{subdir[0][0]}_mra.fasta', 
		      cfg)\
             .posttask(lambda: log_task_progress(cfg, 'assemble_merged', completed=True))\
             .jobs_limit(4)
             #.posttask(lambda: clean_assembly_dir('mra_assembly'))
                    
       
       
            
        # #################################
        #
        #   Q u a l i t y   c o n t r o l
        #
        # #############################


        qc_raw_reads_task = \
            p.transform(qc_fastqs,
                        link_fastqs_task,
                        formatter('.+/(?P<SAMPLE_ID>[^/]+)\.fastq\.gz$'), 
                        os.path.join(cfg.runs_scratch_dir,'qc','read_qc/')+'{SAMPLE_ID[0]}_fastqc.html',
                        cfg,
                        name='qc_raw_reads')\
             .follows(mkdir(os.path.join(cfg.runs_scratch_dir,'qc','read_qc')))


        qc_raw_trimmed_reads_task = \
            p.transform(qc_fastqs, 
                        trim_reads_task,
                        formatter('.+/(?P<SAMPLE_ID>[^/]+)\.fq\.gz$', '.+/(?P<SAMPLE_ID>[^/]+)\.fq\.gz$', None, None), 
                        [os.path.join(cfg.runs_scratch_dir,'qc','read_qc')+'/{SAMPLE_ID[0]}_fastqc.html',
                         os.path.join(cfg.runs_scratch_dir,'qc','read_qc')+'/{SAMPLE_ID[1]}_fastqc.html'],
                        cfg,
                        name='qc_raw_trimmed_reads_task')\
             .follows(mkdir(os.path.join(cfg.runs_scratch_dir,'qc','read_qc')))


        qc_merged_reads_task = \
            p.transform(qc_fastqs,
                        trim_merged_reads_task,
                        formatter('.+/(?P<SAMPLE_ID>[^/]+)\.fq\.gz$'), 
                        os.path.join(cfg.runs_scratch_dir,'qc','read_qc')+'/{SAMPLE_ID[0]}_fastqc.html',
                        cfg,
                        name='qc_merged_reads')\
             .follows(mkdir(os.path.join(cfg.runs_scratch_dir,'qc','read_qc')))
                

        qc_notmerged_pairs_task = \
            p.transform(qc_fastqs,
                        trim_notmerged_pairs,
                        formatter('.+/(?P<SAMPLE_ID>[^/]+)\.fq\.gz$', '.+/(?P<SAMPLE_ID>[^/]+)\.fq\.gz$', None, None), 
                        [os.path.join(cfg.runs_scratch_dir,'qc','read_qc')+'/{SAMPLE_ID[0]}_fastqc.html',
                         os.path.join(cfg.runs_scratch_dir,'qc','read_qc')+'/{SAMPLE_ID[1]}_fastqc.html'],
                        cfg,
                        name = 'qc_notmerged_pairs')\
             .follows(mkdir(os.path.join(cfg.runs_scratch_dir,'qc','read_qc')))


        qc_reads_task = p.follows(name='qc_reads', task_func=do_nothing)\
                         .follows(qc_raw_reads_task, qc_merged_reads_task, qc_notmerged_pairs_task)\
                         .posttask(lambda: log_task_progress(cfg, 'qc_reads', completed=True))




        return p
   
   
