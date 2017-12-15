

import sys
import os
import glob



#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#  Common functions 

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888


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
def run_cmd(cfg, cmd, args, interpreter_args=None, dockerize=None, run_locally=True,
            cpus=1, mem_per_cpu=1024, walltime='24:00:00', 
            retain_job_scripts = True, job_script_dir = None):
    
    if job_script_dir is None:
        job_script_dir = os.path.join(cfg.runs_scratch_dir, "drmaa")
    
    if dockerize is None: dockerize = cfg.dockerize
    
    full_cmd = ("{docker} "+cmd).format(docker = cfg.docker if dockerize else "",
                                        args=args, 
                                        interpreter_args = interpreter_args if interpreter_args!=None else "")

    stdout, stderr = "", ""
    job_options = "--ntasks=1 \
                   --cpus-per-task={cpus} \
                   --mem-per-cpu={mem} \
                   --time={time} \
                  ".format(cpus=cpus, mem=int(1.2*mem_per_cpu), time=walltime)
    #print full_cmd                   
    try:
        stdout, stderr = run_job(full_cmd.strip(), 
                                 job_other_options=job_options,
                                 run_locally = run_locally, 
                                 retain_job_scripts = retain_job_scripts, job_script_directory = job_script_dir,
                                 logger=cfg.logger, working_directory=os.getcwd(),
                                 drmaa_session = cfg.drmaa_session)
    except error_drmaa_job as err:
        raise Exception("\n".join(map(str, ["Failed to run:", cmd, err, stdout, stderr])))


""" 
Currently not available in dockerized mode. 
Only default job scheduling params of run_command available when executing via SLURM.
"""
def run_piped_command(cfg, *args):
    run_locally=True
    retain_job_scripts = True
    job_script_dir = os.path.join(cfg.runs_scratch_dir, "drmaa")	
    cpus=1
    mem_per_cpu=1024
    walltime="24:00:00"
 
    stdout, stderr = "", ""
    job_options = "--ntasks=1 \
                   --cpus-per-task={cpus} \
                   --mem-per-cpu={mem} \
                   --time={time} \
                  ".format(cpus=cpus, mem=int(1.2*mem_per_cpu), time=walltime)
	
    full_cmd = expand_piped_command(*args)
	
    try:
        stdout, stderr = run_job(full_cmd.strip(), 
                                 job_other_options=job_options,
                                 run_locally = run_locally, 
                                 retain_job_scripts = retain_job_scripts, job_script_directory = job_script_dir,
                                 logger=global_vars.cfg.logger, working_directory=os.getcwd(),
                                 drmaa_session = global_vars.cfg.drmaa_session)
    except error_drmaa_job as err:
        raise Exception("\n".join(map(str, ["Failed to run:", full_cmd, err, stdout, stderr])))
	
def expand_piped_command(cmd, cmd_args, interpreter_args=None, *args):
	expanded_cmd = cmd.format(args=cmd_args, interpreter_args = interpreter_args if interpreter_args!=None else "")
	expanded_cmd += (" | "+expand_piped_command(*args)) if len(args) > 0 else ""
	return expanded_cmd


def log_task_progress(cfg, task_name, completed=True):
    cfg.logger.info('Task [%s] %s.' % (task_name, 'completed' if completed else 'started'))


def produce_fastqc_report(cfg, fastq_file, output_dir=None):
    args = fastq_file
    args += (' -o '+output_dir) if output_dir != None else ''
    run_cmd(cfg, cfg.fastqc, args)


def do_nothing():
    pass

    


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   Pipeline tasks

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

from ruffus import *


#
#
# FASTQ generation and archiving
# 

def bcl2fastq_conversion(run_directory, completed_flag, cfg):
    """ Run bcl2fastq conversion and create fastq files in the run directory"""
    out_dir = os.path.join(cfg.runs_scratch_dir,'fastqs')
    interop_dir = os.path.join(out_dir,'InterOp')

    # r, w, and p specify numbers of threads to be used for each of the concurrent subtasks of the conversion (see bcl2fastq manual) 
    args = "-R {indir} -o {outdir} --interop-dir={interopdir} -r1 -w1 -p2 \
            ".format(indir=run_directory, outdir=out_dir, interopdir=interop_dir)
    if cfg.run_on_bcl_tile != None:
        args += " --tiles %s" % cfg.run_on_bcl_tile
        
    run_cmd(cfg, cfg.bcl2fastq, args, cpus=8, mem_per_cpu=2048)
    


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



# -------------------------------------------------------------------- # 


    #8888888888888888888888888888888888888888888888888888
    #
    #            P r e p r o c e s s i n g
    #
    #8888888888888888888888888888888888888888888888888888


#
# Prepare directory for every sample and link the input fastq files
# Expected format:
#    /path/to/file/[SAMPLE_ID]_S[1-9]\d?_L\d\d\d_R[12]_001.fastq.gz
# SAMPLE_ID can contain all signs except path delimiter, i.e. "\"
#
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
def trim_reads(inputs, outfqs, cfg):
    """ Trim reads """
    args = "PE -phred33 -threads 1 \
            {in1} {in2} {out1} {unpaired1} {out2} {unpaired2} \
            ILLUMINACLIP:{adapter}:2:30:10 \
            LEADING:3 \
            TRAILING:3 \
            SLIDINGWINDOW:4:15 \
            MINLEN:36".format(in1=inputs[0], in2=inputs[1],
                              out1=outfqs[0], out2=outfqs[1],
                              unpaired1=outfqs[2], unpaired2=outfqs[3],
                              adapter=adapters)

#    max_mem = 2048
    run_cmd(cfg, cfg.trimmomatic, args) #interpreter_args="-Xmx"+str(max_mem)+"m", cpus=1, mem_per_cpu=max_mem)



#
# Input FASTQ filenames are expected to have following format:
#    [SAMPLE_ID]_[S_NUM]_[LANE_ID]_[R1|R2]_001.fastq.gz
# In this step, overlapping reads from two FASTQ files matching on [SAMPLE_ID] will be merged together. 
# The output will be written to three FASTQ files
#    [SAMPLE_ID]_merged.fq.gz - containing merged reads
#    [SAMPLE_ID]_R1.fq.gz     - with notmerged R1
#    [SAMPLE_ID]_R2.fq.gz     - with notmerged R2 
# SAMPLE_ID can contain all signs except path delimiter, i.e. "\"
# 
def merge_reads(inputs, outputs, cfg):
	""" Merge overlapping reads """
	
	hist=outputs[0].replace('_merged.fq.gz','.hist')
	args='in1={fq1} in2={fq2} \
	      out={fqm} outu1={u1} outu2={u2} \
	      ihist={hist} adapters={adapters} \
	      threads=1'.format(fq1=inputs[0], fq2=inputs[1],
					fqm=outputs[0], u1=outputs[1], u2=outputs[2],
					hist=hist, adapters=cfg.adapters)
		  
	run_cmd(cfg, cfg.bbmerge, args)
    
    
    
#
# Input FASTQ filenames are expected to have following format:
#    [SAMPLE_ID]_R[12].fq.gz
# In this step, the two FASTQ files with nonoverlapping R1 and R2 reads will be trimmed together. 
# The output will be written to two FASTQ files
#    [SAMPLE_ID]_R1.trimmed.fq.gz
#    [SAMPLE_ID]_R2.trimmed.fq.gz
# SAMPLE_ID can contain all signs except path delimiter, i.e. "\"
#
def trim_notmerged_pairs(inputs, outfqs, cfg):
    """ Trim nonoverlapping reads """
    args = "PE -phred33 -threads 1 \
            {in1} {in2} {out1} {unpaired1} {out2} {unpaired2} \
            ILLUMINACLIP:{adapter}:2:30:10 \
            SLIDINGWINDOW:4:15 MINLEN:36 \
            ".format(in1=inputs[1], in2=inputs[2],
                                       out1=outfqs[0], out2=outfqs[1],
                                       unpaired1=outfqs[2], unpaired2=outfqs[3],
                                       adapter=cfg.adapters)
#    max_mem = 2048
    run_cmd(cfg, cfg.trimmomatic, args) #interpreter_args="-Xmx"+str(max_mem)+"m", cpus=1, mem_per_cpu=max_mem)


#
# Input FASTQ filename is expected to have following format:
#    [SAMPLE_ID]_merged.fq.gz
# In this step, the FASTQ file with merged overlapping reads will be trimmed. 
# The output will be written to:
#    [SAMPLE_ID]_merged.trimmed.fq.gz
# SAMPLE_ID can contain all signs except path delimiter, i.e. "\"
#
def trim_merged_reads(input_fqs, trimmed_fq, cfg):
    """ Trim merged overlapping reads """

    merged_fq=input_fqs[0]
    args = "SE -phred33 -threads 1 \
            {fq_in} {fq_out} ILLUMINACLIP:{adapter}:2:30:10 \
            SLIDINGWINDOW:4:15 MINLEN:36 \
            ".format(fq_in=merged_fq, fq_out=trimmed_fq, adapter=cfg.adapters)
#    max_mem = 2048
    run_cmd(cfg, cfg.trimmomatic, args) #interpreter_args="-Xmx"+str(max_mem)+"m", cpus=1, mem_per_cpu=max_mem)




# -------------------------------------------------------------------- # 


    #8888888888888888888888888888888888888888888888888888
    #
    #       M a p p i n g   a n d   c a l l i n g
    #
    #8888888888888888888888888888888888888888888888888888



def bwa_map_and_sort(output_bam, ref_genome, fq1, fq2=None, read_group=None, threads=1):
	
	bwa_args = "mem -t {threads} {rg} {ref} {fq1} \
	            ".format(threads=threads, 
                        rg="-R '%s'" % read_group if read_group!=None else "", 
                        ref=ref_genome, fq1=fq1)
	if fq2 != None:
		bwa_args += fq2

	samtools_args = "sort -o {out}".format(out=output_bam)

	run_piped_command(cfg, cfg.bwa, bwa_args, None,
	                  cfg.samtools, samtools_args, None)

def merge_bams(out_bam, *in_bams):
	threads = 1
	mem = 4096
	
	args = "merge %s" % out_bam
	for bam in in_bams:
		args += (" "+bam)
		
	run_cmd(cfg, cfg.samtools, args, cpus=threads, mem_per_cpu=int(mem/threads))
	
	
def map_reads(fastq_list, ref_genome, output_bam, read_groups=None):
    
    # If no read groups is provided, we could make up default ones based on fq filenames.
    # This would most likely result in unpaired FQs ending up in different read group ids, and in consequence, being genotyped separately
    if read_groups==None:
        s_ids = [ os.path.basename(s[0][:-len('_R1.fq.gz')] if isinstance(s, tuple) else s[:-len('fq.gz')]) for s in fastq_list ]
        read_groups = [ '@RG\tID:{sid}\tSM:{sid}\tLB:{sid}'.format(sid=s) for s in s_ids]
    
    tmp_bams = [ output_bam+str(i) for i in range(0, len(fastq_list)) ]
    for i in range(0, len(fastq_list)):
		if isinstance(fastq_list[i], tuple):
			bwa_map_and_sort(tmp_bams[i], ref_genome, fastq_list[i][0], fastq_list[i][1], read_groups[i])
		else:
			bwa_map_and_sort(tmp_bams[i], ref_genome, fastq_list[i], read_groups[i])   
    
    merge_bams(output_bam, *tmp_bams)
    
    for f in tmp_bams:
		  os.remove(f)


def map_trimmed_reads(fastqs, bam_file, sample_id, cfg):
    """ Maps trimmed paired and unpaired reads. """
    fq1=fastqs[0]
    fq2=fastqs[1]
    fq1u=fastqs[2]
    fq2u=fastqs[3]

    read_groups = ['@RG\tID:{rgid}\tSM:{rgid}\tLB:{lb}'.format(rgid=sample_id, lb=sample_id),
        '@RG\tID:{rgid}\tSM:{rgid}\tLB:{lb}'.format(rgid=sample_id, lb=sample_id+"_U1"),
        '@RG\tID:{rgid}\tSM:{rgid}\tLB:{lb}'.format(rgid=sample_id, lb=sample_id+"_U2"),]

    map_reads([(fq1,fq2),fq1u, fq2u], cfg.reference, bam_file, read_groups)


#
#
#


def call_variants_freebayes(bams_list, vcf, cfg, bam_list_filename='/tmp/bam_list'):
    
    threads = 1
    mem = 4096
    
    with open(bam_list_filename,'w') as f:
        for bam in bams_list:
            f.write(bam + '\n')
    
    args = args = " -f {ref} -v {vcf} -L {bam_list} \
        ".format(ref=cfg.reference, vcf=vcf, bam_list=bam_list_filename)
            
    run_cmd(cfg, cfg.freebayes, args, cpus=threads, mem_per_cpu=int(mem/threads))
    
    os.remove(bam_list_filename)


def call_variants_on_trimmed(bam, vcf, cfg):
    """ Call variants using freebayes on trimmed (not merged) reads """
    call_variants_freebayes([bam], vcf, cfg, bam+'.lst')


def jointcall_variants_on_trimmed(bams, vcf, cfg):
    """ Call variants using freebayes on trimmed (not merged) reads """
    call_variants_freebayes(bams, vcf, cfg)


# -------------------------------------------------------------------- # 


    #8888888888888888888888888888888888888888888888888888
    #
    #                  A s s e m b l y 
    #
    #8888888888888888888888888888888888888888888888888888


def clean_trimmed_fastqs(cfg):
    """ Remove the trimmed fastq files. Links to original fastqs are kept """
    for f in glob.glob(os.path.join(cfg.runs_scratch_dir,'*','*.fq.gz')):
        os.remove(f)

def clean_assembly_dir(assembly_name, cfg):
    """ Remove the temporary assembly files """
    import shutil
    for f in glob.glob(os.path.join(cfg.runs_scratch_dir,'*',assembly_name+'_assembly')):
            print 'rm -r '+f
            #shutil.rmtree(f)


def run_spades(cfg, out_dir, fq=None, fq1=None, fq2=None, 
					fq1_single=None, fq2_single=None, 
					threads = 4, mem_gb=8):
						
    args = "-o {out_dir} -m {mem} -t {threads} --careful ".format(out_dir=out_dir, mem=mem_gb, threads=threads)
	
	# PE inputs if provided
    if fq1 != None and fq2 != None: 
        args += " --pe1-1 {fq1} --pe1-2 {fq2}".format(fq1=fq1,fq2=fq2)
		
	# add SE inputs
    i = 1
    for se_input in [fq, fq1_single, fq2_single]:
        if se_input != None:
            args += " --s{index} {se_input}".format(index=i, se_input=se_input)
            i = i + 1
	
    #print args
    run_cmd(cfg, cfg.spades, args, cpus=threads, mem_per_cpu=int(mem_gb*1024/threads))


def spades_assembly(cfg, fasta_files, assembly_name, **args):
    
    out_dir=os.path.join(os.path.dirname(scaffolds_file), assembly_name)
    if not os.path.isdir(out_dir):
		os.mkdir(out_dir)
        
    run_spades(cfg, out_dir, **args)
    
    import shutil
    shutil.copy(os.path.join(out_dir,'scaffolds.fasta'), fasta_files[0])
    shutil.copy(os.path.join(out_dir,'contigs.fasta'), fasta_files[1])
    #shutil.rmtree(out_dir)


#
# In this step, a pair of trimmed FASTQ files coming from trim_reads are used for assembly.
# The output will be written to SAMPLE_ID directory:
#    [SAMPLE_ID]/tra_assembly
#    [SAMPLE_ID]/[SAMPLE_ID]_tra.fasta
#
def assemble_trimmed(fastqs, scaffolds, cfg):
    fastqs=fastqs[0]   
    spades_assembly(cfg, scaffolds, 'tra_assembly', 
        fq1=fastqs[0], fq2=fastqs[1], 
        fq1_single=fastqs[2], fq2_single=fastqs[3], 
        threads = 4, mem_gb=8)


#
# In this step, the 5 FASTQ files coming from trim_merged_reads and trim_notmerged_pairs are used for assembly.
# The output will be written to SAMPLE_ID directory:
#    [SAMPLE_ID]/mra_assembly
#    [SAMPLE_ID]/[SAMPLE_ID]_mra.fasta
#
def assemble_merged(fastqs, fastas, cfg):
    fqm=fastqs[0]
    fq1=fastqs[1][0]
    fq2=fastqs[1][1]
    fq1u=fastqs[1][2]
    # fq2u is typicaly low quality

    spades_assembly(cfg, fastas, 'mra_assembly', 
        fq=fqm, fq1=fq1, fq2=fq2, 
        fq1_single=fq1u, 
        threads = 4, mem_gb=8)


    
# -------------------------------------------------------------------- #    


    #8888888888888888888888888888888888888888888888888888
    #
    #          Q u a l i t y   c o n t r o l
    #
    #8888888888888888888888888888888888888888888888888888


def qc_fastqs(input_fastqs, reports, cfg):
    """ Produces len(reports) reports for    """
    if not isinstance(input_fastqs, list):
        input_fastqs = [input_fastqs]
    if not isinstance(reports, list):
        reports = [reports]
    if len(input_fastqs)<len(reports):
        raise Exception('Lengths of input FASTQs and report files do not match:\n%s\n%s' % (input_fastqs, reports))
    
    for i in range(0, len(reports)):
        produce_fastqc_report(cfg, input_fastqs[i], os.path.dirname(reports[i]))



def qc_tr_assemblies(scaffolds, report_dir, cfg):
    args = ("-o %s " % report_dir) + " ".join(scaffolds)
    run_cmd(cfg, cfg.quast, args)


def qc_mr_assemblies(scaffolds, report_dir, cfg):
    args = ("-o %s " % report_dir) + " ".join(scaffolds)
    run_cmd(cfg, cfg.quast, args)



# -------------------------------------------------------------------- #


    #8888888888888888888888888888888888888888888888888888
    #
    #     A r c h i v i n g   a n d   c l e a n u p
    #
    #8888888888888888888888888888888888888888888888888888


import shutil

def archive_files(files, archived_paths):
    if not isinstance(files, list):
        files = [files]
    if not isinstance(archived_paths, list):
        archived_paths = [archived_paths]
        
    for i in range(0, len(files)):
        shutil.copyfile(files[i], archived_paths[i])
    
def archive_dir(directory, archived_dir):
    if os.path.exists(archived_dir):
        shutil.rmtree(archived_dir)
    shutil.copytree(directory, archived_dir)  
    #run_cmd(cfg, "cp -r %s %s" % (qc_dir, archived_qc_dir), "", dockerize=False, run_locally=True)

def archive_pipeline(pipeline_dir, archived_paths, cfg):
    run_cmd(cfg, "cp -r %s %s" % (pipeline_dir, archived_paths[0]), "", dockerize=False, run_locally=True)
    shutil.copyfile(pipeline_dir+'.config', archived_paths[1])


def cleanup_files():
    pass
#    run_cmd("rm -rf {dir}/*/mra_assembly \
#            {dir}/*/*merged*.fq.gz \
#            ".format(dir=global_vars.cfg.runs_scratch_dir), "", run_locally=True)





