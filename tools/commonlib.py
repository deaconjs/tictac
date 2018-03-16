import sys
import os
import time
import glob
from datetime import datetime
import imp

homebase = os.getenv('TICTAC_HOME', '.')

# checks and returns input args supplied by user
def check_input(pipeline=None):
    if not pipeline:
        print("Error: no pipeline given to check_input")
        sys.exit()
    if pipeline == 'submitjobs':
        maxjobs = 5
        if len(sys.argv != 2):
            print("usage python submitjobs_script concurrent_threads (e.g. 5)")
            sys.exit()
        try:
            maxjobs = int(sys.argv[1])
        except ValueError:
            print("supply an integer argument e.g. 5")
            sys.exit()
        if maxjobs > 200:
            print("change the code to run %s (>200) threads. setting maxjobs to 200"%(maxjobs))
            maxjobs = 200
        return maxjobs

    if pipeline == 'intervals':
        normalbam, tumorbam, benchmark, genomefile, maxjobs, testorrun = None, None, None, None, 5, 'test'
        if len(sys.argv) != 8:
            print("usage: python run_pipeline.py intervals normalbam tumorbam benchmark genomefile maxjobs testorrun")
            sys.exit()
        # bam files first
        if not os.path.isfile(sys.argv[2]):
            print("supply normal bam or bam.gz file as first argument")
            sys.exit()
        else:
            normalbam = sys.argv[2]
            if not normalbam.endswith('.bam') and not normalbam.endswith('bam.gz'):
                print("Warning: normal bam (arg 1) does not end with .bam or .bam.gz")
                time.sleep(1)
            if 'normal' not in normalbam:
                print("Warning: normal bam (arg 1) does not contain the substring 'normal'")
                time.sleep(1)
            if 'tumor' in normalbam:
                print("Error: Check input arguments, normal bam (arg 1) contains the substring 'tumor'.")
                sys.exit()
        if not os.path.isfile(sys.argv[3]):
            print("supply tumor bam or bam.gz file as second argument")
            sys.exit()
        else:
            tumorbam = sys.argv[3]
            if not tumorbam.endswith('.bam') and not tumorbam.endswith('bam.gz'):
                print("Warning: tumor bam (arg 2) does not end with .bam or .bam.gz")
                time.sleep(1)
            if 'tumor' not in tumorbam:
                print("Warning: tumor bam (arg 2) does not contain the substring 'tumor'")
                time.sleep(1)
            if 'normal' in tumorbam:
                print("Error: Check input arguments, tumor bam (arg 2) contains the substring 'normal'.")
                sys.exit()
        # check benchmark name
        if not os.path.isdir("%s/benchmarks/%s"%(homebase, sys.argv[4])):
            print("supply benchmark name (from tictac/benchmarks e.g. is3) as third argument")
            sys.exit()
        else:
            benchmark = sys.argv[4]
        # check genome file. this can probably be generated silently from the bams' reference .fai if available
        #   instead of requiring user input
        if not os.path.isfile(sys.argv[5]):
            print("supply .genome file (reference .fai columns 1&2) location as fourth argument")
            sys.exit()
        else:
            genomefile = sys.argv[5]
        # check maxjobs for integer and hard limit on number of threads to keep up
        try:
            maxjobs = int(sys.argv[6])
        except ValueError:
            print("supply an integer argument e.g. 5 as fifth argument")
            sys.exit()
        if maxjobs > 200:
            print("change the code to run more than 200 threads.")
            sys.exit()
        # last arg is required declaration of whether this is a test run or full run
        # test simply prints out what its going to do without committing any changes
        testorrun = sys.argv[7]
        if testorrun != "test" and testorrun != "run":
            print("supply 'test' or 'run' as sixth argument")
            sys.exit()
        else:
            if testorrun == 'test':
                print("test mode: splitbam.py will generate SGE submission scripts but not submit them to qsub.")
            else:
                print("run mode: splitbam.py will generate SGE submission scripts and submit to qsub.")
        return normalbam, tumorbam, benchmark, genomefile, maxjobs, testorrun


def list_targets(grep=None, quick=False):
    # list targets available                                                                                                   
    bglb, cglb = [], []
    if grep:
        bglb = glob.glob('%s/benchmarks/*%s*'%(homebase, grep))
        cglb = glob.glob('%s/callers/*%s*'%(homebase, grep))
    else:
        bglb = glob.glob('%s/benchmarks/*'%(homebase))
        cglb = glob.glob('%s/callers/*'%(homebase))

    all_targets_tmp = bglb + cglb
    all_targets = sorted(all_targets_tmp, key=lambda file:os.path.getctime(file))
    all_targets.reverse()
    all_targets_abbrv = []
    for b in bglb + cglb:
        all_targets_abbrv.append(b.split('/')[-1])
    print(time.strftime("%Y-%m-%d %H:%M"))
    print("sh\t|bam\tvcf\tmpil\tmerged\t|r\tqw\t|pipeline")
    print("------------------------------------------------------------------")
    for target in all_targets:
        t = target.split('/')[-1]
        pipeline_dirs_tmp = glob.glob('%s/*'%(target))
        pipeline_dirs_tmp.sort()
        pipeline_dirs = sorted(pipeline_dirs_tmp, key=lambda file:os.path.getctime(file))
        pipeline_dirs.reverse()
        if len(pipeline_dirs) > 0:
            for pd in pipeline_dirs:
                if os.path.isdir(pd):
                    shglob  = glob.glob('%s/work/*.sh'%(pd))
                    pileglob = glob.glob('%s/results/*.mpileup'%(pd))
                    vcfglob = glob.glob('%s/results/*.vcf*'%(pd))
                    if t.startswith('somaticseq'):
                        vcfglob = glob.glob('%s/results/*.vcf/*/*.vcf'%(pd))
                    bamglob = glob.glob('%s/results/*.bam*'%(pd))
                    mergedglob = []
                    if os.path.isdir('%s/final'%(pd)):
                        mergedglob = glob.glob('%s/final/merged*.vcf*'%(pd))
                    p = pd.split('/')[-1]
                    pipelinebase = _get_pipelinebase(t, p)
                    sys.path.append(pipelinebase)
                    import pipeline_commands
                    imp.reload(pipeline_commands)
                    del sys.path[-1]
                    abbv = pipeline_commands.get_abbreviation()
                    running = 0
                    waiting = 0
                    if len(abbv) > 0:
                        cmd = "qstat | grep %s | wc -l > ./scratch/qstatgrepcnt%s.tmp"%(abbv, abbv)
                        os.system(cmd)
                        cmd = "qstat | grep %s > ./scratch/qstatgrep%s.tmp"%(abbv, abbv)
                        os.system(cmd)
                        time.sleep(0.02)
                        f = open('./scratch/qstatgrepcnt%s.tmp'%(abbv), 'r')
                        qstatlines = f.readlines()
                        total = int(qstatlines[0].rstrip('\r\n'))
                        f = open('./scratch/qstatgrep%s.tmp'%(abbv), 'r')
                        qstatlines = f.readlines()
                        time.sleep(0.02)
                        for ql in qstatlines:
                            if 'qw' in ql:
                                waiting += 1
                        running = total - waiting
                        time.sleep(0.02)
                    print("%s\t|%s\t%s\t%s\t%s\t|%s\t%s\t|%s"%(len(shglob), len(bamglob), len(vcfglob), len(pileglob), len(mergedglob), running, waiting, pipelinebase.split('/tictac/')[1]))
        else:
            print("0\t|0\t0\t0\t0\t|0\t0\t|%s"%(t))


# returns number of lines in qstat containing the given caller abbreviation
def qstat_cnt(caller):
    os.system('qstat > qstat_out')
    time.sleep(0.1)
    lines = open('qstat_out', 'r').readlines()
    mycnt = 0
    for line in lines:
        if caller in line:
            mycnt += 1
    return mycnt


# checks input arguments against directory/file structure and returns pipeline base
def _get_pipelinebase(target, pipeline):
    refbase="%s/reference"%(homebase)
    # check input
    # find the target string in benchmarks and callers
    benchglob = glob.glob('%s/benchmarks/*'%(homebase))
    callerglob = glob.glob('%s/callers/*'%(homebase))
    targetbase = None
    for f in benchglob:
        if f.endswith('/%s'%(target)):
            targetbase = f
    for f in callerglob:
        if f.endswith('/%s'%(target)):
            targetbase = f
    if not targetbase:
        print("target %s not found in benchmarks/ or callers/. please check input"%(target))
        sys.exit()
    targetglob = glob.glob('%s/*'%(targetbase))
    pipelinebase = '%s/%s'%(targetbase, pipeline)
    if not os.path.isdir(pipelinebase):
        print("pipeline directory %s not found. To run pipeline %s create that directory and populate with pipeline_commands.py and results/, logs/, and work/ directories."%(pipelinebase, pipeline))
        sys.exit()
    if not os.path.isfile('%s/pipeline_commands.py'%(pipelinebase)):
        print("pipeline_commands.py not found in pipeline directory %s. add a pipeline_commands.py file there, with get_run_command() and get_extension() and get_abbreviation() function declared. see other pipelines for examples."%(pipelinebase))
        sys.exit()
    run_command = '%s/pipeline_commands.py'%(pipelinebase)
    run_command_lines = open(run_command, 'r').readlines()
    # more check input
    for line in run_command_lines:
        if "def get_run_command" in line:
            break
    else:
        print("no 'def get_run_command' in %s/pipeline_commands.py. add a get_run_command that returns the run commands to add to a qsub script. See other pipelines for examples."%(pipelinebase))
        sys.exit()
    return pipelinebase


def train_caller(target, pipeline, maxjobs=10):
    print("training caller %s using benchmark %s"%(target, pipeline))
    pipelinebase = _get_pipelinebase(target, pipeline)
    if not os.path.isdir('%s/training'%(pipelinebase)):
        os.system('mkdir %s/training'%(pipelinebase))
        os.system('mkdir %s/training/work'%(pipelinebase))
        os.system('mkdir %s/training/results'%(pipelinebase))
        os.system('mkdir %s/training/final'%(pipelinebase))

    if not os.path.isfile('%s/training_commands.py'%(pipelinebase)):
        print("no training_commands.py file found in %s. Creating from template."%(pipelinebase))
        training_templates_dir = 'tools/std_templates/training_templates'
        if os.path.isfile('%s/%s.py'%(training_templates_dir, target)):
            os.system('cp %s/%s.py %s/training_commands.py'%(training_templates_dir, target, pipelinebase))
        else:
            print("no standard template %s.py available in %s. please create one from the other examples there to train %s on %s"%(target, training_templates_dir, target, pipeline))
            sys.exit()
        print("Customize the training_commands.py file in %s (generated from %s/%s.py)\n"%(pipelinebase, training_templates_dir, target))
        sys.exit()

    print("caller %s ready for training"%(target))

    build_submission_scripts(target, pipeline, True)
    run_submission_scripts(target, pipeline, maxjobs, 'run', True)
    merge_training_results(target, pipeline)

    
def merge_training_results(target, pipeline):
    pipelinebase = _get_pipelinebase(target, pipeline)
    sys.path.append(pipelinebase)
    import training_commands
    imp.reload(training_commands)
    del sys.path[-1]
    training_commands.merge_results(target, pipeline)

# Searches for target directory in benchmarks and callers, then search for
# subdiretory by pipeline name, which should have a pipeline_commands.py file with
# a get_run_command.py function in it.
#
# e.g. build_submission_scripts('is3', 'mpileups', '.mpileup')
# e.g. build_submission_scripts('is3', 'splitbam', '.bam.gz')
# e.g. build_submission_scripts('mutect', 'is3', '.vcf')
#
def build_submission_scripts(target, pipeline, istraining=False):
    pipelinebase = _get_pipelinebase(target, pipeline)
    print("building submission scripts for %s %s at %s"%(target, pipeline, pipelinebase))
    # first delete previous submission scripts
    # and import or reload training/pipeline commands
    sys.path.append(pipelinebase)
    import pipeline_commands
    imp.reload(pipeline_commands)
    if istraining:
        os.system('rm %s/training/work/*'%(pipelinebase))
        import training_commands
        imp.reload(training_commands)
    else:
        os.system('rm %s/work/*'%(pipelinebase))
    del sys.path[-1]

    extension = pipeline_commands.get_extension()
    sampletypes = ['tumor', 'normal']   # ['combined']
    sampletypes = pipeline_commands.get_sampletypes()
    refbase="%s/reference"%(homebase)
    # sort out the regions, extract and sort indices       
    glb = glob.glob("%s/bedfiles/*.bed"%(refbase))
    print(len(glb), " bed regions")
    indices = {}
    for f in glb:
        indices[int(f.split('/')[-1].split('.')[0])] = f
    keys = list(indices.keys())
    keys.sort()
    # build the submission scripts
    for x, index in enumerate(keys):
        lines = open(indices[index], 'r').readlines()
        fname = "%s.bed"%(x+1)
        bedname = '%s/bedfiles/%s'%(refbase, fname)
        for sampletype in sampletypes:
            if istraining:
                outlines = training_commands.get_training_command(target, pipeline, index, sampletype)
            else:
                outlines = pipeline_commands.get_run_command(index, sampletype)
            if not istraining:
                subfilename = "%s/work/%s.%s.%s.%s.sh"%(pipelinebase, sampletype, target, pipeline, index)
            else:
                subfilename = "%s/training/work/%s.%s.%s.%s.sh"%(pipelinebase, sampletype, target, pipeline, index)
            subfile = open(subfilename, 'w')
            subfile.writelines(outlines)
            subfile.close()


def submit(script, pipelinebase):
    # log time and script contents                                                                     
    outlines = []
    outlines.append("\nqueue submission %s"%(time.strftime("%Y-%m-%d %H:%M")))
    for line in open(script, 'r').readlines():
        outlines.append('\t%s'%(line))
    f = open('./log.log', 'a')
    f.writelines(outlines)
    f.close()
    os.system('chmod 777 %s'%(script))
    os.system('qsub -e %s/logs -o %s/logs %s'%(pipelinebase, pipelinebase, script))


def run_submission_scripts(target, pipeline, maxjobs, mode, istraining=False):
    pipelinebase = _get_pipelinebase(target, pipeline)
    sys.path.append(pipelinebase)
    import pipeline_commands
    imp.reload(pipeline_commands)
    abbreviation = pipeline_commands.get_abbreviation()
    if istraining:
        import training_commands
        imp.reload(training_commands)
        abbreviation = training_commands.get_abbreviation()
        
    del sys.path[-1]
    extension = pipeline_commands.get_extension()
    sampletypes = pipeline_commands.get_sampletypes()
    glb = None
    if not istraining:
        glb = glob.glob("%s/work/*.sh"%(pipelinebase))
    else:
        glb = glob.glob("%s/training/work/*.sh"%(pipelinebase))
    glb.sort()
    print("running %s threads, max %s at a time"%(len(glb), maxjobs))
    sys.stdout.flush()
    maxjobsname = None
    if not istraining:
        maxjobsname = './maxjobs.%s.%s.tmp'%(target, pipeline)
    else:
        maxjobsname = './maxjobs.%s.%s.training.tmp'%(target, pipeline)
    f = open(maxjobsname, 'w')
    f.write('%s\n'%(maxjobs))
    f.close()
    time.sleep(0.1)

    msg = "tick tock "   # each cycle 5 min
    cycle_time = 30      # seconds
    msgcnt = 0
    print("current qstat shows %s jobs with abbreviation %s"%(qstat_cnt(abbreviation), abbreviation))
    submitted_jobs_cnt = 0
    submitted_jobs_indices = []
    if maxjobs:
        for script in glb:
            # read from maxjobs file; user can adjust the maximum number of jobs by editing this file          
            while qstat_cnt(abbreviation) >= int(maxjobs):
                # see if user has adjusted maxjobs
                new_maxjobs = int(open(maxjobsname,'r').readlines()[0].rstrip('\r\n'))
                time.sleep(0.1)
                if maxjobs != new_maxjobs:
                    print("\nuser adjusted maxjobs. Old %s, new %s."%(maxjobs, new_maxjobs))
                    maxjobs = new_maxjobs
                    # hard limit. Adjust this to allow over 200 concurrent threads at a time                    
                    if maxjobs > 200:
                        print("to set maxjobs over 200 edit the code. setting maxjobs to 200.")
                        maxjobs = 200
                    sys.stdout.flush()
                    if new_maxjobs > maxjobs: # ok to go right to submission                                   
                        break
                # write timer message                                                                          
                sys.stdout.write('%s'%(msg[msgcnt%len(msg)]))
                sys.stdout.flush()
                msgcnt += 1
                time.sleep(cycle_time)
            tokens = script.split('/')
            index = tokens[-1].split('.')[3]
            sampletype = tokens[-1].split('.')[0]
            expout, ss_expout = None, None
            if not istraining:
                expout = '%s/results/%s.%s.%s.%s.%s'%(pipelinebase, sampletype, target, pipeline, index, extension)
            else:
                expout = '%s/training/results/%s.%s.%s.%s.%s'%(pipelinebase, sampletype, target, pipeline, index, extension)
            # for somaticseq
            if not istraining:
                ss_expout = '%s/results/%s.%s'%(pipelinebase, index, extension)
            else:
                ss_expout = '%s/training/results/%s.%s'%(pipelinebase, index, extension)
            if os.path.isfile(expout) or os.path.isfile('%s.snp'%(expout)) or os.path.isfile('%s.indel'%(expout)) or os.path.isdir(ss_expout):
                shortexpout = expout.split(pipelinebase)[1]
                if os.path.isdir(ss_expout):
                    shortexpout = ss_expout.split(pipelinebase)[1]
                print("already have %s, skip submission for %s"%(shortexpout, script.split('/')[-1]))
                sys.stdout.flush()
            else:
                print("%s submitting script %s"%(mode, script))
                sys.stdout.flush()
                submitted_jobs_cnt += 1
                submitted_jobs_indices.append(index)
                if mode == 'run':
                    try:
                        t = int(maxjobs)
                    except TypeError:
                        print("non-integer argument (%s) given for maxjobs. check input and resubmit"%(maxjobs))
                    else:
                        maxjobs = t
                    submit(script, pipelinebase)
                    msgcnt = 0
                    time.sleep(2)

    os.system('rm ./maxjobs.%s.%s.tmp'%(target, pipeline))
    print("%s jobs submitted for mode %s, target %s, pipeline %s, indices %s"%(submitted_jobs_cnt, mode, target, pipeline, ','.join(submitted_jobs_indices)))

# if mode is inspect just print if mode is 'fix' then fix selected
def check_submission_results(target, pipeline, mode='inspect', fixtype=None):
    pipelinebase = _get_pipelinebase(target, pipeline)
    sys.path.append(pipelinebase)
    import pipeline_commands
    imp.reload(pipeline_commands)
    del sys.path[-1]
    abbreviation = pipeline_commands.get_abbreviation()
    extension = pipeline_commands.get_extension()
    print("checking in %s/results/*.%s"%(pipelinebase, extension))
    snpglb = glob.glob('%s/results/*.%s.snp.vcf'%(pipelinebase, extension))
    indelglb = glob.glob('%s/results/*.%s.indel.vcf'%(pipelinebase, extension))
    tumorglb = glob.glob('%s/results/*tumor.*.%s'%(pipelinebase, extension))
    normalglb = glob.glob('%s/results/*normal.*.%s'%(pipelinebase, extension))
    vcfglb = glob.glob('%s/results/*.%s'%(pipelinebase, extension))
    # for somaticseq
    aredirs = False
    for v in vcfglb:
        if v.endswith('.%s'%(extension)) and os.path.isdir(v):
            aredirs = True
    if aredirs:
        print("looking in %s/results/*.%s/*SNV*.%s"%(pipelinebase, extension, extension))
        snpglb = glob.glob('%s/results/*.%s/*SNV*.%s'%(pipelinebase, extension, extension)) + glob.glob('%s/results/*.%s/*/*SNV*.%s'%(pipelinebase, extension, extension))
        indelglb = glob.glob('%s/results/*.%s/*INDEL*.%s'%(pipelinebase, extension, extension))
        print("found %s snp vcfs, %s indel vcfs"%(len(snpglb), len(indelglb)))
        vcfglb = []

    snpglb.sort()     # all of the *.vcf.snp files from results unless SS then *SNV*.vcf
    indelglb.sort()   # all of the *.vcf.indel files from results unless SS then *INDEL*.vcf
    tumorglb.sort()   # all of the *tumor.*.vcf files from results
    normalglb.sort()  # all of the *normal.*.vcf files from results
    vcfglb.sort()     # all of the *.vcf files from results unless SS then empty

    sampletype = 'other'    
    glb = vcfglb
    if len(tumorglb) > 0 and len(normalglb) > 0:
        sampletype = 'tumornormal'
        glb = tumorglb + normalglb
    if len(snpglb) > 0 or len(indelglb) > 0:
        sampletype = 'snpindel'
        glb = snpglb #+ indelglb

    print("running sampletype %s with %s combined samples"%(sampletype, len(glb)))

    # glb now contains the target vcf file list
    log_files_dir = '%s/logs'%(pipelinebase)
    log_glb = glob.glob('%s/*'%(log_files_dir))

    index_has_trained = {}
    for vcffile in glb:
        ##print "## vcffile %s"%(vcffile)
        vcffiledir  = '/'.join(vcffile.split('/')[:-1])
        # is somaticseq?
        index = None
        if aredirs:
            index = ""
            if vcffiledir.endswith('.vcf'):
                index = int(vcffiledir.split('/')[-1].split('.')[0])
            else:
                index = int(vcffiledir.split('/')[-2].split('.')[0])
            if index not in list(index_has_trained.keys()):
                index_has_trained[index] = False
            #print vcffile
            if vcffile == '%s/Trained.sSNV.vcf'%(vcffiledir):
                index_has_trained[index] = True
                #print 'index %s val %s'%(index, index_has_trained[index])
            else:
                #print 'index %s val %s'%(index, index_has_trained[index])
                continue
        else:
            index = int(vcffile.split('/')[-1].split('.')[3])

        #print "querying %s/%s%s.e*"%(log_files_dir, abbreviation, index)
        thisvcfelogs = glob.glob('%s/%s%s.e*'%(log_files_dir, abbreviation, index)) + glob.glob('%s/%s%scombined.e*'%(log_files_dir, abbreviation, index))
        thisvcfologs = glob.glob('%s/%s%s.o*'%(log_files_dir, abbreviation, index)) + glob.glob('%s/%s%scombined.o*'%(log_files_dir, abbreviation, index))
        if len(thisvcfelogs) == 0 and len(thisvcfologs) == 0:
            thisvcfelogs = glob.glob('%s/%s%snormal.e*'%(log_files_dir, abbreviation, index)) + glob.glob('%s/%s%stumor.e*'%(log_files_dir, abbreviation, index))
            thisvcfologs = glob.glob('%s/%s%snormal.o*'%(log_files_dir, abbreviation, index)) + glob.glob('%s/%s%stumor.e*'%(log_files_dir, abbreviation, index))
            
        if len(thisvcfelogs) > 1:
            #print "%s in thisvcfelogs, %s in thisvcfologs"%(len(thisvcfelogs), len(thisvcfologs))
            #print "%s, %s"%(thisvcfelogs[0], thisvcfologs[0])
            thisvcfelogs.sort(key=os.path.getmtime)
            thisvcfologs.sort(key=os.path.getmtime)
            thisvcfelogs.reverse()
            thisvcfologs.reverse()
            thisvcfelogs = [thisvcfelogs[0]]
            thisvcfologs = [thisvcfologs[0]]

        if len(thisvcfologs) > 0 and len(thisvcfelogs) > 0:
            #print '%s e logs %s o logs'%(thisvcfologs, thisvcfelogs)
            cmd = 'echo %s; tail -n 2 %s'%(thisvcfelogs[0].split('/')[-1], thisvcfelogs[0])
            os.system(cmd)
            cmd = 'echo %s; tail -n 2 %s'%(thisvcfologs[0].split('/')[-1], thisvcfologs[0])
            os.system(cmd)
        else:
            print("!*!*!*! WARNING:::: no log files found for %s/%s%s"%(log_files_dir, abbreviation, index))

    if aredirs:
        keys = list(index_has_trained.keys())
        keys.sort()
        kmin, kmax = keys[0], keys[-1]
        missing_trained_count = 0        
        for key in range(kmin, kmax+1):
            if key not in index_has_trained or index_has_trained[key] == False:
                print("no trained results vcf found for index %s"%(key))
                missing_trained_count += 1
        if missing_trained_count == 0:
            print("Training results found for all vcfs")

    # report any files more than n SD above or % of median
    sizes = {}
    for vcffile in glb:
        sizes[vcffile] = os.path.getsize(vcffile)
    sum = 0
    for vcffile in glb:
        sum += sizes[vcffile]
    avg = int(sum/len(glb))
    sum = 0
    for vcffile in glb:
        sum += (sizes[vcffile]-avg)**2
    sd = int(sum/(len(glb)-1))**0.5
    justsizes = []
    for key in list(sizes.keys()):
        justsizes.append(sizes[key])
    justsizes.sort()
    med = justsizes[int(len(justsizes)/2)]

    fs_zero = []
    kb=1024
    print('median %s kb'%(int(med/kb)))
    print('avg %s kb'%(int(avg/kb)))
    print('sd  %s kb'%(int(sd/kb)))
    smallest = {}
    for vcffile in glb:
        vcffile_local = vcffile.split(homebase)[1]
        if not list(smallest.keys()):
            smallest[vcffile] = sizes[vcffile]
        if sizes[vcffile] < smallest[list(smallest.keys())[0]]:
            smallest[vcffile] = sizes[vcffile]
        if sizes[vcffile] == 0:
            fs_zero.append(vcffile)
        if sizes[vcffile] > avg+(5*sd):
            print("size %s kb > 5xSD above avg for %s"%(int(sizes[vcffile]/kb), vcffile_local))
        if sizes[vcffile] < 0.01*med:
            if sizes[vcffile] < kb:
                print("size %s bytes < 1%% group median for %s"%(sizes[vcffile], vcffile_local))
            else:
                print("size %s kb < 1%% group median for %s"%(int(sizes[vcffile]/kb), vcffile_local))
    print("%s files are size zero"%(len(fs_zero)))
    print("smallest bam %s %s"%(list(smallest.keys())[0], smallest[list(smallest.keys())[0]]))

    missing_i = []
    missing_t = []
    missing_n = []
    missing_s = []
    missing_i = []
    if sampletype == 'tumornormal':
        # report any missing files                                                                              
        t_indices = []
        n_indices = []
        for vcffile in glb:
            tokens = vcffile.split('/')[-1].split('.')
            if 'tumor' in tokens:
                t_indices.append(int(tokens[3]))
            elif 'normal' in tokens:
                n_indices.append(int(tokens[3]))
        t_indices.sort()
        N_indices.sort()

        t_max, t_min = max(t_indices), min(t_indices)
        n_max, n_min = max(n_indices), min(n_indices)

        for x, index in enumerate(range(t_min, t_max+1)):
            if index not in t_indices:
                missing_t.append(index)
        for x, index in enumerate(range(n_min, n_max+1)):
            if index not in n_indices:
                missing_n.append(index)

        if len(missing_t) == 0:
            print("tumor indices appear to be complete")
        else:
            missing_t.sort()
            print("missing tumor indices:")
            print(missing_t)

        if len(missing_n) == 0:
            print("normal indices appear to be complete")
        else:
            missing_t.sort()
            print("missing normal indices:")
            print(missing_n)
    elif sampletype == 'snpindel':
        # report any missing files                                                                              
        s_indices = []
        i_indices = []
        for vcffile in glb:
            tokens = vcffile.split('/')[-1].split('.')
            if tokens[-1] == 'snp':
                s_indices.append(int(tokens[3]))
            elif 'SNV' in tokens[1]:
                s_indices.append(int(vcffile.split('/')[-2].split('.')[0]))
            elif tokens[-1] == 'indel':
                i_indices.append(int(tokens[3]))
            elif 'INDEL' in tokens[1]:
                i_indices.append(int(vcffile.split('/')[-2].split('.')[0]))
        s_indices.sort()
        i_indices.sort()

        # look for snp index ranges
        if len(s_indices) > 0:
            s_max, s_min = max(s_indices), min(s_indices)
            for x, index in enumerate(range(s_min, s_max+1)):
                if index not in s_indices:
                    missing_s.append(index)
            if len(missing_s) == 0:
                print("snp indices appear to be complete")
            else:
                missing_s.sort()
                print("missing snp indices:")
                print(missing_s)

        # look for indel index ranges
        if len(i_indices) > 0:
            i_max, i_min = max(i_indices), min(i_indices)
            for x, index in enumerate(range(i_min, i_max+1)):
                if index not in i_indices:
                    missing_i.append(index)
            if len(missing_i) == 0:
                print("indel indices appear to be complete")
            else:
                missing_i.sort()
                print("missing indel indices:")
                print(missing_i)
        
    elif sampletype == 'other':
        indices = []
        for vcffile in glb:
            tokens = vcffile.split('/')[-1].split('.')
            indices.append(int(tokens[3]))
        indices.sort()
        imax, imin = max(indices), min(indices)
        for x, index in enumerate(range(imin, imax+1)):
            if index not in indices:
                missing_i.append(index)
        if len(missing_i) == 0:
            print("indices appear to be complete")
        else:
            missing_i.sort()
            print("missing indices:")
            print(missing_i)


    if mode == 'fix' and fixtype == 'missing_indices' and (len(missing_i) != 0 or len(missing_t) != 0 or len(missing_n) != 0 or len(missing_s) != 0):
        print('missing indices indicates bam generation failed. Fix any issues with manual testing on qsub and direct command line calls, then Re-issue the run commmand to attempt to fill the gaps, and check again.')
        sys.exit()

    if len(fs_zero):
        print("files with size zero:")
        for fsz in fs_zero:
            print(fsz)
        if mode == 'fix' and fixtype == 'zero_file_size':
            for fsz in fs_zero:
                print('removing size zero bam file %s'%(fsz))
                os.system('rm %s'%(fsz))
            print("removed %s size zero bam files. Now re-issue the run command to fill those gaps, then check again."%(len(fs_zero)))
    else:
        print("no files found with size zero")


def create_target(target, caller_or_benchmark):
    target_dir = '%s/%ss/%s'%(homebase, caller_or_benchmark, target)
    if not os.path.isdir(target_dir):
        os.system('mkdir %s'%(target_dir))
    time.sleep(0.1)
    """
    """

def create_benchmark(target):
    create_target(target, 'benchmark')

def create_caller(target):
    create_target(target, 'caller')

def extract_pipeline(target, template, callerorbenchmark):
    target_dir = '%s/%ss/%s'%(homebase, callerorbenchmark, target)
    pipeline_dir = "%s/%s"%(target_dir, template)
    templates_dir = '%s/tools/std_templates/%s_templates'%(homebase, callerorbenchmark)
    pipeline_templates = glob.glob('%s/*.tar'%(templates_dir))

    for ct in pipeline_templates:
        if ct.endswith('.tar') and template == ct.split('/')[-1].split('.tar')[0]:
           #if not os.path.isdir(pipeline_dir):
           #    os.system('mkdir %s'%(pipeline_dir))
           cmd = 'cp %s/%s.tar %s'%(templates_dir, template, target_dir)
           print(cmd)
           os.system(cmd)
           os.system('tar -xvf %s/%s.tar -C %s'%(target_dir, template, target_dir))
           os.system('rm %s/%s.tar'%(target_dir, template))
           break
    else:
        print("no tarballs in %s match template %s"%(templates_dir, template))

    print("now edit the file %s/pipeline_commands.py to populate abbreviation, extension, and command line details"%(pipeline_dir))
    sys.exit()


def _run_merge(target, pipeline, glb, merged_vcf_out):
    try:
        f=open(merged_vcf_out, 'w')
        f.writelines(["\n"])
        f.close()
    except:
        print("run_merge.py cannot write to output file %s"%(merged_vcf_out))
        sys.exit()

    # build vcf.list file
    vcffiles = {}
    if len(glb) == 0:
        print("no vcf files provided, please check input")
        sys.exit()

    outlines = []
    for vcff in glb:
        if not vcff.endswith('.gz'):
            outlines.append('%s.gz\n'%(vcff))
        else:
            outlines.append('%s\n'%(vcff))
    vcf_list_file = open('./scratch/%s.%s.vcf.list'%(target, pipeline), 'w')
    vcf_list_file.writelines(outlines)
    vcf_list_file.close()
    time.sleep(0.1)

    os.system('. /etc/profile.d/modules.sh')
    os.system('module load bgzip')
    os.system('module load tabix')

    # compress with bgzip
    print("compressing with bgzip")
    for vcff in glb:
        #print "looking for %s.gz"%(vcff)
        if not vcff.endswith('.gz') and not os.path.isfile("%s.gz"%(vcff)):
            bcmd = 'bgzip -c %s > %s.gz'%(vcff, vcff)
            print(bcmd)
            os.system(bcmd)

    # index with tabix
    print("indexing with tabix")
    for vcff in glb:
        if not os.path.isfile("%s.tbi"%(vcff)) and not os.path.isfile("%s.gz.tbi"%(vcff)):
            tcmd = ""
            if not vcff.endswith('.gz'):
                tcmd = 'tabix -p vcf %s.gz'%(vcff)
            else:
                tcmd = 'tabix -p vcf %s'%(vcff)
            print(tcmd)
            os.system(tcmd)

    bcfcmd = "bcftools concat "
    bcfcmd += "--file-list ./scratch/%s.%s.vcf.list "%(target, pipeline)
    bcfcmd += "--output %s "%(merged_vcf_out)
    bcfcmd += "--output-type v "
    bcfcmd += "--threads 1 "
    bcfcmd += "-a "
    print("running %s"%(bcfcmd))
    os.system('module load bcftools python/2.7.8; %s'%(bcfcmd))


def merge_vcfs(target, pipeline):
    pipelinebase = _get_pipelinebase(target, pipeline)

    # grab files by combined.target.pipeline.###.vcf
    vcf_dir  = "%s/results"%(pipelinebase)

    # glob has some weird regex that requires iteratively adding digits per order of magnitude
    # automatically checks for .gz versions
    def ranged_glob(strlist1, strlist2, magnitude):
        holder = []
        for i in range(1, magnitude):
            holder += glob.glob('.'.join(strlist1 + '[0-9]'*i + strlist2))
        for i in range(1, magnitude):
            holder += glob.glob('.'.join(strlist1 + '[0-9]'*i + strlist2 + ['.gz']))
        return sorted(holder, key=lambda x: int(x.split('.')[-len(strlist2)]))

    aredirs = False      # are the vcfs directories?

    caller = None
    known_callers = ['varscan', 'lofreq', 'muse', 'somaticseq', 'strelka', 'sniper', 'vardict', 'mutect', 'mutect2', 'jointsnvmix']
    for c in known_callers:
        if target.startswith(c):
            caller = c

    vcfglb = []  # original vcfs
    svcfglb = [] # sorted vcfs

    # look for varscan first 
    vcfglb = ranged_glob(['%s/combined'%(vcf_dir), target, pipeline], ['vcf', 'snp', 'vcf'], 5) # look for unsorted
    svcfglb = ranged_glob(['%s/combined'%(vcf_dir), target, pipeline], ['vcf', 'snp', 'sorted', 'vcf'], 5) # and sorted
    if len(vcfglb) >  0 or len(svcfglb) > 0:
        caller = 'varscan'

    if not caller:
        vcfglb = ranged_glob(['%s/combined'%(vcf_dir), target, pipeline], ['vcfsomatic_final_minus-dbsnp', 'snvs', 'vcf'], 5)
        svcfglb = ranged_glob(['%s/combined'%(vcf_dir), target, pipeline], ['vcfsomatic_final_minus-dbsnp', 'snvs', 'sorted', 'vcf'], 5)
        if len(vcfglb) > 0 or len(svcfglb) > 0:
            caller = 'lofreq'

    if not caller: 
        # cover normal case including mutect, muse, sniper, etc
        vcfglb = ranged_glob(['%s/combined'%(vcf_dir), target, pipeline], ['vcf'], 5)
        svcfglb = ranged_glob(['%s/combined'%(vcf_dir), target, pipeline], ['sorted', 'vcf'], 5)
        caller = 'general'    

    if len(vcfglb) == 0 and len(svcfglb) == 0:
        tvcfglb = glob.glob('%s/*.vcf'%(vcf_dir))
        vcfglb = sorted(tvcfglb, key=lambda x:int(x.split('/')[-1].split('.')[0]))
        caller = 'unknown'

    print("caller %s covers %s files, %s sorted'%(caller, len(vcfglb), len(svcfglb)))


    # handle strelka and somaticseq first
    aredirs = False
    for t in vcfglb:
        if t.endswith('vcf') and os.path.isdir(t):
            aredirs = True  

    if aredirs and (target.startswith('SS') or caller == 'somaticseq'):
        # first determine whether trained or not
        ss_trained = False
        snpglb = glob.glob('%s/*.vcf/*SNV*.vcf'%(vcf_dir)) + glob.glob('%s/*.vcf/*/*SNV*.vcf'%(vcf_dir))
        indelglb = glob.glob('%s/*.vcf/*INDEL*.vcf'%(vcf_dir)) + glob.glob('%s/*.vcf/*/*INDEL*.vcf'%(vcf_dir))
        vcfglb = snpglb + indelglb
        for tglb in vcfglb:
            if "Trained" in tglb:
                ss_trained = True
        # and get the files from trained or untrained file handles
        if ss_trained:
            snpglb = glob.glob('%s/*.vcf/Trained*SNV*.vcf'%(vcf_dir)) + glob.glob('%s/*.vcf/*/Trained*SNV*.vcf'%(vcf_dir))
            indelglb = glob.glob('%s/*.vcf/Trained*INDEL*.vcf'%(vcf_dir)) + glob.glob('%s/*.vcf/*/Trained*INDEL*.vcf'%(vcf_dir))
        else:
            snpglb = glob.glob('%s/*.vcf/Untrained*SNV*.vcf'%(vcf_dir))
            indelglb = glob.glob('%s/*.vcf/Untrained*INDEL*.vcf'%(vcf_dir))
        if len(snpglb) > 0:
            merged_vcf_out = "%s/final/merged.%s.%s.snp.vcf"%(pipelinebase, target, pipeline)
            _run_merge(target, pipeline, snpglb, merged_vcf_out)
        if len(indelglb) > 0:
            merged_vcf_out = "%s/final/merged.%s.%s.indel.vcf"%(pipelinebase, target, pipeline)
            _run_merge(target, pipeline, indelglb, merged_vcf_out)
        sys.exit()

    if aredirs and caller == 'strelka':
        snpglb = glob.glob('%s/*.vcf/results/variants/somatic.snvs.vcf.gz'%(vcf_dir))
        indelglb = glob.glob('%s/*.vcf/results/variants/somatic.indels.vcf.gz'%(vcf_dir))
        if len(snpglb) > 0:
            merged_vcf_out = "%s/final/merged.%s.%s.snp.vcf"%(pipelinebase, target, pipeline)
            _run_merge(target, pipeline, snpglb, merged_vcf_out)
        if len(indelglb) > 0:
            merged_vcf_out = "%s/final/merged.%s.%s.indel.vcf"%(pipelinebase, target, pipeline)
            _run_merge(target, pipeline, indelglb, merged_vcf_out)
        sys.exit()
        
    # next sort any unsorted vcfs
    orig_indices   = {}
    sorted_indices = {}
    index_index = -2
    if caller == 'varscan':
        index_index = -4
    if caller == 'lofreq':
        index_index = -5
    for v in vcfglb:
        orig_indices[int(v.split('.')[index_index])] = v
    for sv in svcfglb:
        sorted_indices[int(sv.split('.')[index_index-1])] = sv
    orig_keys = list(orig_indices.keys()).sort()
    sorted_keys = list(sorted_indices.keys()).sort()

    seqdict = '%s/reference/genomefiles/b37.dict'%(homebase)
    if caller in ['mutect', 'mutect2']:
        seqdict = '%s/reference/genomefiles/b37_decoy.dict'%(homebase)

    for orig_key in range(min(orig_keys), max(orig_keys)+1):
        if orig_key not in orig_keys:
            print("no unsorted version of index %s found, skipping"%(orig_key))
            continue
        if orig_key not in sorted_keys:
            # sort the vcf file
            # start by fixing unruly vcfs
            if caller == 'muse':
                outlines = []
                lines = open(orig_indices[orig_key], 'r').readlines()
                seen_NC = False
                for line in lines:
                    if "NC_007605" in line:
                        if not seen_NC:
                            outlines.append(line)
                            outlines.append("##contig=<ID=hs37d5,length=35477943>\n")
                            seen_NC = True
                    elif "hs37d5" in line:
                        # already added it above so no need to re-add the original
                        pass
                    else:
                        outlines.append(line)
                outf = open("%s.tmp"%(orig_indices[orig_key]), 'w')
                outf.writelines(outlines)
                outf.close()
                time.sleep(0.01)
                os.system('mv %s.tmp %s'%(orig_indices[orig_key], orig_indices[orig_key]))

            elif caller == 'varscan':
                outlines = []
                lines = open(orig_indices[orig_key], 'r').readlines()
                found_nonstandard = False
                in_calls = False
                for line in lines:
                    if line.startswith('#CHROM'):
                        in_calls = True
                        continue
                    if not in_calls:
                        continue
                    if line.split()[3] in ['A', 'T', 'G', 'C', 'N', 'a', 't', 'c', 'g', 'n']:    # case insensitive
                        outlines.append(line)
                    else:
                        found_nonstandard = True
                        print("removing nonstandard REF line %s"%(line))
                if found_nonstandard:
                    outf = open("%s.tmp"%(orig_indices[orig_key]), 'w')
                    outf.writelines(outlines)
                    outf.close()
                    time.sleep(0.05)
                    os.system('mv %s.tmp %s'%(orig_indices[orig_key], orig_indices[orig_key]))
            elif caller == 'jointsnvmix' or target.startswith('jsm'):
                tmpname = '%s.tmp'%('.'.join(orig_indices[orig_key].split('.')[:-1]))
                os.system('mv %s %s'%(orig_indices[orig_key], tmpname))
                time.sleep(0.05)
                os.system('%s/reference/somaticseq/JSM2VCF.sh %s > %s'%(homebase, tmpname, orig_indices[orig_key]))


            # rebuild the sorted filename
            tokens = orig_indices[orig_key].split('.')
            sorted_fname = '%s.sorted.vcf'%('.'.join(tokens[:-1]))
            cmd = 'java -jar %s/tools/picard/picard.jar SortVcf I=%s O=%s SEQUENCE_DICTIONARY=%s'%(homebase,  orig_indices[orig_key], sorted_fname, seqdict)
            print(cmd)
            os.system(cmd)

    mergevcfglb   = []
    mergeindelglb = []
    mergesnpglb   = []

    for tv in vcfglb:
        sorted_tv = '%s.sorted.vcf'%('.'.join(tv.split('.')[:-1]))
        if os.path.isfile(sorted_tv) and os.path.getsize(sorted_tv) > 0:
            mergevcfglb.append(sorted_tv)
        else:
            print("file %s does not exist or is size zero"%(tv))
        if 'indel' in sorted_tv.split('.'):
            mergeindelglb.append(sorted_tv)
        if 'snp' in sorted_tv.split('.'):
            mergesnpglb.append(sorted_tv)

    # check for output directory and contents
    if not os.path.isdir(vcf_dir):
        print("vcf_dir %s must point to a directory with decompressed .vcf files to merge"%(vcf_dir))
        sys.exit()

    if len(mergevcfglb) == 0 and len(mergeindelglb) == 0 and len(mergesnpglb) == 0:
        print("vcf_dir %s must contain some decompressed .vcf and/or .snp.vcf and .indel.vcf files to merge"%(vcf_dir))
        sys.exit()

    # and merge
    if len(mergevcfglb) > 0:
        merged_vcf_out = "%s/final/merged.%s.%s.vcf"%(pipelinebase, target, pipeline)
        _run_merge(target, pipeline, mergevcfglb, merged_vcf_out)

    if len(mergesnpglb) > 0:
        merged_vcf_out = "%s/final/merged.%s.%s.snp.vcf"%(pipelinebase, target, pipeline)
        _run_merge(target, pipeline, mergesnpglb, merged_vcf_out)

    if len(mergeindelglb) > 0:
        merged_vcf_out = "%s/final/merged.%s.%s.indel.vcf"%(pipelinebase, target, pipeline)
        _run_merge(target, pipeline, mergeindelglb, merged_vcf_out)


def compare_vcfs_to_truth(target, pipeline, truthvcf, variant_class=None):
    pipelinebase = _get_pipelinebase(target, pipeline)
    # the merged vcfs
    callervcf = "%s/final/merged.%s.%s.vcf"%(pipelinebase, target, pipeline)
    snpcallervcf = "%s/final/merged.%s.%s.snp.vcf"%(pipelinebase, target, pipeline)
    indelcallervcf = "%s/final/merged.%s.%s.indel.vcf"%(pipelinebase, target, pipeline)

    vcfevaluatorresultsbasename, snpevaluatorresultsbasename, indelevaluatorresultsbasename = "", "", ""
    if os.path.isfile(callervcf):
        vcfevaluatorresultsbasename = "%s/final/merged.%s.%s.evaluator"%(pipelinebase, target, pipeline)
    if os.path.isfile(snpcallervcf):
        snpevaluatorresultsbasename = "%s/final/merged.%s.%s.snp.evaluator"%(pipelinebase, target, pipeline)
    if os.path.isfile(indelcallervcf):
        indelevaluatorresultsbasename = "%s/final/merged.%s.%s.indel.evaluator"%(pipelinebase, target, pipeline)

    eval_combined_SNV_file = "%s.SNV"%(vcfevaluatorresultsbasename)
    eval_combined_INDEL_file = "%s.INDEL"%(vcfevaluatorresultsbasename)
    eval_snp_SNV_file = "%s.SNV"%(snpevaluatorresultsbasename)
    eval_indel_INDEL_file = "%s.INDEL"%(indelevaluatorresultsbasename)

    # print the files if they exist
    previous_results_exist = False
    for r in [eval_combined_SNV_file, eval_snp_SNV_file]:
        if os.path.isfile(r):
            print("\nSNVs")
            os.system('cat %s'%(r))
            previous_results_exist = True
    for r in [eval_combined_INDEL_file, eval_indel_INDEL_file]:
        if os.path.isfile(r):
            print("\nINDELs")
            os.system('cat %s'%(r))
            previous_results_exist = True

    # if previous results do exist, ask whether to overwrite
    if previous_results_exist:
        good_input = False
        while not good_input:
            rawinput = input("\noverwrite previous results (y/n)? ")
            if rawinput not in ['y', 'n']:
                print("%s given. please input y or n"%(rawinput))
                continue
            if rawinput == 'n':
                sys.exit()
            elif rawinput == 'y':
                good_input = True

    variant_classes = []
    if variant_class and variant_class not in variant_classes:
        print("fatal error: variant class %s not found for compare mode"%(variant_class))
        sys.exit()
    if not variant_class:
        variant_classes = ['SNV', 'INDEL'] # , 'SV']
    else:
        variant_classes = [variant_class]

    print("looking for file %s"%(callervcf))
    if os.path.isfile(callervcf):
        for vc in variant_classes:
            cmd = "module load python/2.7.8; python ./tools/dream_evaluator.py %s %s %s > %s.%s"%(callervcf, truthvcf, vc, vcfevaluatorresultsbasename, vc)
            print("\nrunning %s"%(cmd))
            os.system(cmd)
            time.sleep(0.1)
            os.system('cat %s.%s'%(vcfevaluatorresultsbasename, vc))

    if os.path.isfile(snpcallervcf) and 'SNV' in variant_classes:
        cmd = "module load python/2.7.8; python ./tools/dream_evaluator.py %s %s SNV > %s"%(snpcallervcf, truthvcf, eval_snp_SNV_file)
        print("\nrunning %s"%(cmd))
        os.system(cmd)
        time.sleep(0.1)
        os.system('cat %s'%(eval_snp_SNV_file))

    if os.path.isfile(indelcallervcf) and 'INDEL' in variant_classes:
        cmd = "module load python/2.7.8; python ./tools/dream_evaluator.py %s %s INDEL > %s"%(indelcallervcf, truthvcf, eval_indel_INDEL_file)
        print("\nrunning %s"%(cmd))
        os.system(cmd)
        time.sleep(0.1)
        os.system('cat %s'%(eval_indel_INDEL_file))


def caller_report(masked_or_unmasked, mut_type, greps=[]):
    caller_dirs = []
    if len(greps) > 0:
        # get all the callers that have a grep in their name or in one of their pipeline names
        for grep in greps:
            caller_dirs = caller_dirs + glob.glob('%s/callers/*%s*'%(homebase, grep))
            subdirs = glob.glob('%s/callers/*/*%s*'%(homebase, grep))
            for subdir in subdirs:
                cdir = '/'.join(subdir.split('/')[:-1])
                if cdir not in caller_dirs:
                    caller_dirs.append(cdir)
    else:
        caller_dirs = glob.glob('%s/callers/*'%(homebase))

    # get the latest timestamp from vcfs for each pipeline
    pipeline_times = {}
    times = []
    for caller in caller_dirs:
        if caller not in pipeline_times:
            pipeline_times[caller] = {}
        pipelines = glob.glob('%s/*'%(caller))
        latest_caller_time = os.path.getctime(caller)
        for pipeline in pipelines:
            pipeline_times[caller][pipeline] = os.path.getctime("%s/results"%(pipeline))
            times.append(pipeline_times[caller][pipeline])

    times.sort()
    times.reverse()
    sorted_caller_pipeline_pairs = []
    for time in times:
        for caller in list(pipeline_times.keys()):
            for pipeline in list(pipeline_times[caller].keys()): 
                if pipeline_times[caller][pipeline] == time:
                    sorted_caller_pipeline_pairs.append([caller, pipeline])

    pipeline_dirs = []
    print("\ntictac caller report %s %s\n"%(masked_or_unmasked, mut_type))
    print("tp#\tfp#\tsubrecs\tsubmask\ttrurecs\ttrumask\tsens\tspec\tb_accu\ttotal\tlastrun\tcaller\tpipeline")
    print("------------------------------------------------------------------------")
    for caller_pipeline_pair in sorted_caller_pipeline_pairs:
        caller_dir, pipeline_dir = caller_pipeline_pair[0], caller_pipeline_pair[1]
        caller_glb = glob.glob('%s/*'%(caller_dir))
        caller = caller_dir.split('/')[-1].rstrip('\r\n')
        if os.path.isdir(pipeline_dir):
            pipeline = pipeline_dir.split('/')[-1].rstrip('\r\n')
            if len(greps) > 0:
                for grep in greps:
                    if grep in pipeline or grep in caller:
                        break
                else:
                    continue
            final_glb = glob.glob('%s/final/*'%(pipeline_dir))
            for final_file in final_glb:
                if final_file.endswith('.%s'%(mut_type)):
                    lines = open(final_file, 'r').readlines()
                    if len(lines) < 12:
                        print('check %s/%s/final for evaluator.py results'%(caller, pipeline))
                        continue
                    masked_counts_line = lines[3].rstrip('\r\n')
                    masked_stats_line = lines[4].rstrip('\r\n')
                    masked_tots_line = lines[5].rstrip('\r\n')
                    unmasked_counts_line = lines[9].rstrip('\r\n')
                    unmasked_stats_line = lines[10].rstrip('\r\n')
                    unmasked_tots_line = lines[11].rstrip('\r\n')
                    (mtpc, mfpc, msre, msma, mtre, mtma) = masked_counts_line.rstrip('\r\n').split(' ')
                    (a, b, c, mspe, mbal) = masked_stats_line.rstrip('\r\n').split(',')
                    (d, e, f, msen) = c.split(' ')
                    mnum = masked_tots_line.split(' ')[-1]
                    (utpc, ufpc, usre, usma, utre, utma) = unmasked_counts_line.rstrip('\r\n').split(' ')
                    (a, b, c, uspe, ubal) = unmasked_stats_line.rstrip('\r\n').split(',')
                    (d, e, f, usen) = c.split(' ')
                    unum = unmasked_tots_line.split(' ')[-1]
                    date = datetime.fromtimestamp(pipeline_times[caller_dir][pipeline_dir]).strftime('%Y-%m-%d')
                    if masked_or_unmasked == 'masked':
                        print("%s\t%s\t%s\t%s\t%s\t%s\t%4.3f\t%4.3f\t%4.3f\t%s\t%s\t%s\t%s"%(mtpc, mfpc, msre, msma, mtre, mtma, float(msen), float(mspe), float(mbal), mnum, date, caller, pipeline))
                    elif masked_or_unmasked == 'unmasked':
                        print("%s\t%s\t%s\t%s\t%s\t%s\t%4.3f\t%4.3f\t%4.3f\t%s\t%s\t%s\t%s"%(utpc, ufpc, usre, usma, utre, utma, float(usen), float(uspe), float(ubal), unum, date, caller, pipeline))

"""
masked:
tpcount, fpcount, subrecs, submasked, trurecs, trumasked:
3845 195 4040 2683 3979 353
sensitivity, specificity, balanced accuracy: 0.966323196783,0.951732673267,0.959027935025
number of unmasked mutations in submission: 4040

unmasked:
tpcount, fpcount, subrecs, submasked, trurecs, trumasked:
4178 217 4395 1839 4322 10
sensitivity, specificity, balanced accuracy: 0.966682091624,0.950625711035,0.95865390133
number of unmasked mutations in submission: 4395
"""
def watch_queue(grep=None):
    monitoring = True
    grep = ""
    print_list = True
    while monitoring:
        if print_list:
            list_targets(grep)
        print_list = True
        grep = ""
        available_actions = {'':'continue', 'q':'quit monitoring', 'r':'report', 'substring':'grep'} #, 'p':'pause pipeline', 't':'terminate pipeline', 'm':'adjust maxjobs'}
        good_input = False
        action = ''
        rawwatchinput = ""
        while not good_input:
            query = ''
            for aakey in list(available_actions.keys()):
                query += '(%s)%s, '%(aakey, available_actions[aakey])
            rawwatchinput = input('which action %s : '%(query))
            if True: #rawwatchinput in available_actions.keys():
                action = rawwatchinput
                print("taking action %s"%(available_actions[action]))
                good_input = True
            #else:
            #    print "(%s) not among %s. check input"%(rawwatchinput, ','.join(available_actions.keys()))
        if action == 'q':
            break
        elif action == '':
            continue
        elif action == 'r':
            caller_report('masked', 'SNV')
            print_list = False
            action = ''
            continue 
        else:
            grep = rawwatchinput
        #print "taking action %s"%(available_actions[action])


