        
from ruffus import *
from pipeline.tasks import *
    
class ZaoPipelineFactory:    
    
    def __init__(self):
        pass
    
    
    def get_fastq_generation_subpipeline(self, cfg, pipeline):
        """ Generate FASTQ, archive them and link into folders """
        
        bcl2fastq_conversion_task = pipeline\
            .transform(bcl2fastq_conversion, 
                        cfg.run_folder, formatter(), 
                        os.path.join(cfg.runs_scratch_dir,'fastqs','completed'), 
                        cfg,
                        name='bcl2fastq_conversion')\
            .follows(mkdir(cfg.runs_scratch_dir))\
            .follows(mkdir(os.path.join(cfg.runs_scratch_dir,'fastqs')))\
            .posttask(touch_file(os.path.join(cfg.runs_scratch_dir,'fastqs','completed')))\
            .posttask(lambda: log_task_progress(cfg, 'bcl2fastq_conversion', completed=True))


        archive_fastqs_task = pipeline\
            .transform(archive_fastqs, 
                       bcl2fastq_conversion_task,
                       formatter(".+/(?P<RUN_ID>[^/]+)/fastqs/completed"), 
                       str(cfg.fastq_archive)+"/{RUN_ID[0]}/fastq",
                       name = 'archive_fastqs')\
            .active_if(cfg.fastq_archive != None)\
            .posttask(lambda: log_task_progress(cfg, 'archive_fastqs', completed=True))
    

        link_fastqs_task = pipeline\
            .transform(link_fastqs,
                       cfg.input_fastqs,
                       formatter('(?P<PATH>.+)/(?P<SAMPLE_ID>[^/]+)_S[1-9]\d?_L\d\d\d_R[12]_001\.fastq\.gz$'),
                       cfg.runs_scratch_dir+'/{SAMPLE_ID[0]}/{basename[0]}{ext[0]}',
                       name = 'link_fastqs')\
            .follows(archive_fastqs_task)\
            .jobs_limit(1)\
            .posttask(lambda: log_task_progress(cfg, 'link_fastqs', completed=True))


        qc_raw_reads_task = pipeline\
            .transform(qc_fastqs,
                       link_fastqs_task,
                       formatter('.+/(?P<SAMPLE_ID>[^/]+)\.fastq\.gz$'), 
                       os.path.join(cfg.runs_scratch_dir,'qc','read_qc/')+'{SAMPLE_ID[0]}_fastqc.html',
                       cfg,
                       name='qc_raw_reads')\
            .follows(mkdir(os.path.join(cfg.runs_scratch_dir,'qc','read_qc')))


        archive_raw_reads_qc_task = pipeline\
            .transform(archive_files,
                   qc_raw_reads_task,
                   formatter(),
                   os.path.join(cfg.results_archive, cfg.run_id, 'qc', 'read_qc', '{basename[0]}{ext[0]}'),
                   name = 'archive_raw_read_qc')\
            .follows(mkdir(os.path.join(cfg.results_archive, cfg.run_id, 'qc', 'read_qc')))\
            .active_if(cfg.results_archive != None)

        return pipeline
    

    def get_merge_and_trim_subpipeline(self, cfg, pipeline, inputs):
        """ Merge overlapping reads, and trim the merged and not merged reads. INPUTS: output_from('link_fastqs')"""
        
        merge_reads_task = \
            pipeline.collate(merge_reads,
                             inputs,
                             regex(r'(.+)/([^/]+)_S[1-9]\d?_(L\d\d\d)_R[12]_001\.fastq\.gz$'), 
                             [r'\1/\2_merged.fq.gz', r'\1/\2_notmerged_R1.fq.gz', r'\1/\2_notmerged_R2.fq.gz'],
                             cfg,
                             name = 'merge_reads')\
                    .posttask(lambda: log_task_progress(cfg, 'merge_reads', completed=True))


        trim_notmerged_pairs_task = \
            pipeline.transform(trim_notmerged_pairs, 
                               merge_reads_task,
                               formatter(None, '.+/(?P<PREFIX>[^/]+)\.fq\.gz$', '.+/(?P<PREFIX>[^/]+)\.fq\.gz$'), 
                               ['{path[1]}/{PREFIX[1]}.trimmed.fq.gz', '{path[2]}/{PREFIX[2]}.trimmed.fq.gz',
                                '{path[1]}/{PREFIX[1]}.unpaired.fq.gz', '{path[2]}/{PREFIX[2]}.unpaired.fq.gz'],
                                cfg,
                               name = 'trim_notmerged_pairs')\
                    .posttask(lambda: log_task_progress(cfg, 'trim_notmerged_pairs', completed=True))
                                                
        
        
        trim_merged_reads_task = \
            pipeline.transform(trim_merged_reads, 
                               merge_reads_task,
                               suffix('_merged.fq.gz'), 
                               '_merged.trimmed.fq.gz',
                               cfg,
                               name = 'trim_merged_reads')\
                    .posttask(lambda: log_task_progress(cfg, 'trim_merged_reads', completed=True))
                    

        qc_merged_reads_task = pipeline\
            .transform(qc_fastqs,
                       trim_merged_reads_task,
                       formatter('.+/(?P<SAMPLE_ID>[^/]+)\.fq\.gz$'), 
                       os.path.join(cfg.runs_scratch_dir,'qc','read_qc')+'/{SAMPLE_ID[0]}_fastqc.html',
                       cfg,
                       name='qc_merged_reads')\
            .follows(mkdir(os.path.join(cfg.runs_scratch_dir,'qc','read_qc')))
                

        qc_notmerged_pairs_task = pipeline\
            .transform(qc_fastqs,
                       trim_notmerged_pairs_task,
                       formatter('.+/(?P<SAMPLE_ID>[^/]+)\.fq\.gz$', '.+/(?P<SAMPLE_ID>[^/]+)\.fq\.gz$', None, None), 
                       [os.path.join(cfg.runs_scratch_dir,'qc','read_qc')+'/{SAMPLE_ID[0]}_fastqc.html',
                        os.path.join(cfg.runs_scratch_dir,'qc','read_qc')+'/{SAMPLE_ID[1]}_fastqc.html'],
                       cfg,
                       name = 'qc_notmerged_pairs')\
            .follows(mkdir(os.path.join(cfg.runs_scratch_dir,'qc','read_qc')))
        
        
        archive_merged_read_qc_task = pipeline\
            .transform(archive_files,
                       qc_merged_reads_task,
                       formatter(),
                       os.path.join(cfg.results_archive, cfg.run_id, 'qc', 'read_qc', '{basename[0]}{ext[0]}'),
                       name = 'archive_merged_read_qc')\
            .follows(mkdir(os.path.join(cfg.results_archive, cfg.run_id, 'qc', 'read_qc')))\
            .active_if(cfg.results_archive != None)
        
        archive_notmerged_pairs_qc_task = pipeline\
            .transform(archive_files,
                       qc_notmerged_pairs_task,
                       formatter(),
                       [os.path.join(cfg.results_archive, cfg.run_id, 'qc', 'read_qc', '{basename[0]}{ext[0]}'),
                        os.path.join(cfg.results_archive, cfg.run_id, 'qc', 'read_qc', '{basename[1]}{ext[1]}')],
                       name = 'archive_notmerged_pairs_qc')\
            .follows(mkdir(os.path.join(cfg.results_archive, cfg.run_id, 'qc', 'read_qc')))\
            .active_if(cfg.results_archive != None)
        
        return pipeline


    def get_trim_subpipeline(self, cfg, pipeline, inputs):
        """ Trim read pairs. INPUTS: output_from('link_fastqs')"""
        
        if isinstance(inputs, basestring):
            inputs = output_from(inputs)

        trim_reads_task = \
            pipeline.collate(trim_reads, 
                             inputs,
                             regex(r'(.+)/([^/]+)_S[1-9]\d?_(L\d\d\d)_R[12]_001\.fastq\.gz$'),
                             [r'\1/\2_\3_R1.fq.gz', r'\1/\2_\3_R2.fq.gz', 
                              r'\1/\2_\3_R1_unpaired.fq.gz', r'\1/\2_\3_R2_unpaired.fq.gz'],
                             cfg,
                             name = 'trim_reads')\
                    .posttask(lambda: log_task_progress(cfg, 'trim_reads', completed=True))
        
        qc_trimmed_reads_task = pipeline\
            .transform(qc_fastqs, 
                       trim_reads_task,
                       formatter('.+/(?P<SAMPLE_ID>[^/]+)\.fq\.gz$', '.+/(?P<SAMPLE_ID>[^/]+)\.fq\.gz$', None, None), 
                       [os.path.join(cfg.runs_scratch_dir,'qc','read_qc')+'/{SAMPLE_ID[0]}_fastqc.html',
                        os.path.join(cfg.runs_scratch_dir,'qc','read_qc')+'/{SAMPLE_ID[1]}_fastqc.html'],
                       cfg,
                       name='qc_trimmed_reads_task')\
            .follows(mkdir(os.path.join(cfg.runs_scratch_dir,'qc','read_qc')))

        return pipeline

    
    def get_mapping_and_jointcalling_subpipeline(self, cfg, pipeline, inputs):
        """ Map and call variants. INPUTS: output_from('trim_reads')"""
        
        if isinstance(inputs, basestring):
            inputs = output_from(inputs)
        
        map_trimmed_reads_task = \
            pipeline.transform(map_trimmed_reads,
                                inputs,
                                formatter("(.+)/(?P<SAMPLE_ID>[^/]+)_L\d\d\d_R[12](_unpaired)?\.fq\.gz$"),
                                "{subpath[0][0]}/{subdir[0][0]}.bam", "{SAMPLE_ID[0]}",
                                cfg,
                                name = 'map_trimmed_reads')\
                    .posttask(lambda: log_task_progress(cfg, 'map_trimmed_reads', completed=True))
   
   
        call_variants_on_trimmed_task = \
            pipeline.transform(call_variants_on_trimmed,
                                map_trimmed_reads_task,
                                suffix(".bam"), 
                                ".fb.vcf",
                                cfg,
                                name = 'call_variants_on_trimmed')\
                    .posttask(lambda: log_task_progress(cfg, 'call_variants_on_trimmed', completed=True))

   
        joincall_variants_on_triummed_task = \
            pipeline.merge(jointcall_variants_on_trimmed,
                            map_trimmed_reads_task, 
                            os.path.join(cfg.runs_scratch_dir, "multisample.fb.vcf"),
                            cfg,
                            name = 'jointcall_variants_on_trimmed')\
                    .posttask(lambda: log_task_progress(cfg, 'jointcall_variants_on_trimmed', completed=True))
   
        return pipeline
        
        
   
    def get_assembly_subpipeline(sefl, cfg, pipeline, inputs):
        """ Assemble merged trimmed reads. INPUTS: [output_from('trim_merged_reads'), output_from('trim_notmerged_pairs')] """
        
        assemble_merged_task = \
            pipeline.collate(assemble_merged,
                             inputs, 
                             formatter(), 
                             ['{subpath[0][0]}/{subdir[0][0]}_mra.fasta', '{subpath[0][0]}/{subdir[0][0]}_mra.fasta.contigs.fasta'], 
                             cfg,
                             name = 'assemble_merged')\
                    .posttask(lambda: log_task_progress(cfg, 'assemble_merged', completed=True))\
                    .jobs_limit(4)
                    #.posttask(lambda: clean_assembly_dir('mra_assembly'))

        archive_fasta_task = pipeline\
            .transform(archive_files,
                       assemble_merged_task,
                       formatter(), 
                       [os.path.join(cfg.results_archive, cfg.run_id, 'fasta', '{basename[0]}{ext[0]}'),
                        os.path.join(cfg.results_archive, cfg.run_id, 'fasta', '{basename[1]}{ext[1]}')],
                       name = 'archive_fasta')\
            .follows(mkdir(os.path.join(cfg.results_archive, cfg.run_id, 'fasta')))\
            .active_if(cfg.results_archive != None)
       
        qc_mr_assemblies_task = pipeline\
            .merge(qc_mr_assemblies,
                   assemble_merged_task,
                   os.path.join(cfg.runs_scratch_dir, 'qc', 'assembly_qc','mr_report'),
                   cfg,
                   name = 'qc_mr_assemblies')\
            .follows(mkdir(os.path.join(cfg.runs_scratch_dir,'qc','assembly_qc')))\
            .posttask(lambda: log_task_progress(cfg, 'qc_mr_assemblies', completed=True))

        archive_assembly_qc_task = pipeline\
            .transform(archive_dir,
                       qc_mr_assemblies_task,
                       formatter(), 
                       os.path.join(cfg.results_archive, cfg.run_id, 'qc','assembly_qc','{basename[0]}'),
                       name = 'archive_assembly_qc')\
            .follows(mkdir(os.path.join(cfg.results_archive, cfg.run_id, 'qc', 'assembly_qc')))\
            .active_if(cfg.results_archive != None)

        return pipeline

    
    
    
    def get_bacterial_assembly_pipeline(self, name, cfg):
        
        # assemble the pipeline 
        p = Pipeline(name=name)
        
        p = self.get_fastq_generation_subpipeline(cfg, p)
        p = self.get_merge_and_trim_subpipeline(cfg, p, p['link_fastqs'])
        p = self.get_assembly_subpipeline(cfg, p, [p['trim_merged_reads'], p['trim_notmerged_pairs']])
        
        archive_pipeline_task = \
            p.files(archive_pipeline,
                    os.path.join(cfg.runs_scratch_dir, 'pipeline'),
                    [os.path.join(cfg.results_archive, cfg.run_id, 'pipeline'), 
                     os.path.join(cfg.results_archive, cfg.run_id, 'pipeline', 'pipeline.config')],
                    cfg)\
             .follows(mkdir(os.path.join(cfg.results_archive, cfg.run_id)))\
             .active_if(cfg.results_archive != None)
       
        # wrap it up
        complete_run_task = \
            p.follows(name='complete_run', task_func=do_nothing)\
             .follows(p["archive_raw_read_qc"])\
             .follows(p["archive_merged_read_qc"])\
             .follows(p["archive_notmerged_pairs_qc"])\
             .follows(p["archive_fasta"])\
             .follows(p["archive_assembly_qc"])\
             .follows(archive_pipeline_task)\
             .posttask(cleanup_files)\
             .posttask(lambda: log_task_progress(cfg, 'complete_run', completed=True))

        return p
   
   
