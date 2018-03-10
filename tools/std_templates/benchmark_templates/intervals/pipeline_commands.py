import sys

def get_extension():
    return "bam.gz"

def get_abbreviation():
    # e.g. return amlinv 

def get_sampletypes():
    return ['normal', 'tumor']

def get_run_command(index, sampletype, queue=None):
    homebase = os.getenv('TICTAC_HOME', '.')
    if not queue:
        myqueue = os.getenv('TICTAC_QUEUE')
    extension = get_extension()
    try:
        abbreviation = get_abbreviation()
    except:
        print("abbreviation has not been set. Have you customized the loaded default pipeline_commands.py file?")
        sys.exit()

    sampletypes = get_sampletypes()

    if not extension or not abbreviation or not sampletypes:
        print "need to add extension, abbreviation, and sampletypes to pipeline_commands.py"
        sys.exit()
    # format input bam files
    bam_dir   =      # e.g. "%s/reference/aml_bams"%(homebase)
    normalbam =      # e.g. "%s/AML31_normal_wgs_na_hiseq_bwa-0.5.9.bam"%(bam_dir)
    tumorbam  =      # e.g. "%s/AML31_primary-tumor_wgs_na_hiseq_bwa-0.5.9.bam"%(bam_dir)

    # format output file
    bench_dir =      # e.g. "%s/benchmarks"%(homebase)
    out_dir   =      # e.g. "%s/aml/intervals/results"%(bench_dir)
    out_file  =      # e.g. "%s/%s.aml.intervals.%s.%s"%(out_dir, sampletype, index, extension) 

    # first two columns from .fai file from reference genome tell bedtools the chromosomal order
    genomefile =     # e.g. "%s/reference/genomefiles/b37.genome"%(homebase)
    # format bed file 
    bed_dir =        # e.g. "%s/reference/bedfiles"%(homebase)
    bed_file =       # e.g. "%s/%s.bed"%(bed_dir, index)

    outlines = ["#!/bin/bash\n", "#$ -S /bin/bash\n", "#$ -N %s%s%s\n"%(abbreviation, index, sampletype), "#$ -q %s\n"%(myqueue), "#$ -cwd\n"] #, ". /etc/profile.d/modules.sh\n"]
    #outlines.append("module load bedtools\n")
    bamtosplit = normalbam
    if sampletype == 'tumor' or sampletype == 'tumour':
        bamtosplit = tumorbam
    cmd = "bedtools intersect -sorted -g %s -b %s -abam %s > %s\n"%(genomefile, bed_file, bamtosplit, out_file)
    outlines.append(cmd)
    return outlines
