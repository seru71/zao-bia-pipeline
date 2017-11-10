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



class SampleTable:
    
    def get_sample_configs:
        # returns a dict with sample_ID : JSON 



if __name__ == '__main__':


    # run bcl2fastq first?

    st = SampleTable()
    cfgs = st.get_sample_configs()

    # group by configs
    
     for each group:
        
       cfg = PipelineConfig(logger)
       cfg.set_runfolder(path)
       cfg.load_from_file(settings_template)
       cfg.load_from_JSON(groups_settings)
       
       import ruffus
       from pipeline import tasks
       
       pipeline_run(cfg.task, jobs=cfg.jobs, .... )
       
       
       
