import sys
sys.path.append('%s/tools'%(os.getenv('TICTAC_HOME', '.'))
import commonlib
import glob
import os
import time



def main(argv):
    # if no arguments, then go into interactive mode. list the available pipeline components 
    #   and offer entry criteria identical to that required on the command line. For each 
    #   command executed, print out the command that would be required at the prompt to 
    #   generate the same behavior. The sequence of commands at the interactive prompt 
    #   should be identical to the arguments in the call.
    #
    # upon submission inform the user on queue monitoring and policy, what happens if the
    #   jobs are interrupted and how to restart, how to check, fix, and merge the final vcfs
    #
    # operational modes: prepare(?), build, run, check, fix, merge, compare
    # information avail: caller command line sets/parameters, statistics on reference or
    #    sample genomes, file QC?

    #home_base = '/home/dsweeney/tools/tictac'
    home_base = os.getenv('TICTAC_HOME', '.')
    interactive_mode = False
    command_line_args = []

    print(argv)
    if len(argv) == 1:
        print("entering interactive mode")
        interactive_mode = True

    # input argument 1: mode
    modes = ['list', 'watch', 'train', 'build', 'run', 'check', 'fix', 'merge', 'compare', 'report']
    if not interactive_mode:
        mode = sys.argv[1]
        if mode not in modes:
            print("mode %s not available. check input"%(mode))
            sys.exit()
    else:
        good_input = False
        while not good_input:
            print("modes %s available"%(', '.join(modes)))
            rawinput = input("which mode? ")
            if rawinput not in modes:
                print("%s not in available modes"%(rawinput))
            else:
                mode = rawinput
                good_input = True

    command_line_args.append(mode)

    if mode == 'list':
        if len(sys.argv) == 2:
            commonlib.list_targets()
            sys.exit()
        elif len(sys.argv) == 3:
            commonlib.list_targets(sys.argv[2])
            sys.exit()

    if mode == 'watch':
        if len(sys.argv) == 2:
            commonlib.watch_queue()
            sys.exit()
        elif len(sys.argv) == 3:
            commonlib.watch_queue(sys.argv[2])
            sys.exit()

    if mode == 'report':
        print(sys.argv)
        # summary report on all evaluator results
        if len(sys.argv) == 2:
            commonlib.caller_report('masked', 'SNV')
            commonlib.caller_report('unmasked', 'SNV')
        elif len(sys.argv) == 3:
            commonlib.caller_report('masked', 'SNV', sys.argv[2:])
            commonlib.caller_report('unmasked', 'SNV', sys.argv[2:])
        sys.exit()

    # check input argument 2: target 
    target = None
    callerorbenchmark = None
    if not interactive_mode:
        target = sys.argv[2]
        target_base = None
    else:
        good_input = False
        while not good_input:
            # list targets available
            # commonlib.list_targets()
            bglb = glob.glob('%s/benchmarks/*'%(home_base))
            cglb = glob.glob('%s/callers/*'%(home_base))
            all_targets = bglb + cglb
            all_targets_abbrv = []
            for b in bglb + cglb:
                all_targets_abbrv.append(b.split('/')[-1])
            print("targets %s available"%(', '.join(all_targets_abbrv)))
            rawtargetinput = ""
            while len(rawtargetinput) < 1:
                rawtargetinput = input("which target? ")
            if len(rawtargetinput.strip()) == 0 or rawtargetinput not in all_targets_abbrv:
                print("%s not in available targets"%(rawtargetinput))
                rawcreateinput = input('target %s not found in benchmarks/ or callers/. add a (b)enchmark, (c)aller, ()none: '%(rawtargetinput))
                if rawcreateinput in ['b', 'c']:
                    if rawcreateinput == 'b':
                        commonlib.create_benchmark(rawtargetinput)
                    elif rawcreateinput == 'c':
                        commonlib.create_caller(rawtargetinput)
                else:
                    print('none selected. hit b or c instead to start a new target')
            else:
                target = rawtargetinput
                good_input = True

    bglb = glob.glob('%s/benchmarks/*'%(home_base))
    cglb = glob.glob('%s/callers/*'%(home_base))
    if '%s/benchmarks/%s'%(home_base, target) in bglb:
        print("target benchmark %s found"%(target))
        callerorbenchmark = 'benchmark'
        target_base = '%s/benchmarks/%s'%(home_base, target)
    elif '%s/callers/%s'%(home_base, target) in cglb:
        print("target caller %s found"%(target))
        callerorbenchmark = 'caller'
        target_base = '%s/callers/%s'%(home_base, target)
    else:
        print('target %s not found in benchmarks/ or callers/.'%(target))
        sys.exit()

    command_line_args.append(target)

    # check input argument 3: pipeline    
    pipeline = None
    if not interactive_mode:
        pipeline = sys.argv[3]
        pipeline_base = None
    else:
        good_input = False
        while not good_input:
            # first find pipelines already available
            glb = glob.glob('%s/*'%(target_base))
            pipeline_names = []
            for g in glb:
                if os.path.isdir(g):
                    pipeline_names.append(g.split('/')[-1])
            # now look through templates for any that can be added
            glb = glob.glob('%s/tools/std_templates/%s_templates/*.tar'%(home_base, callerorbenchmark))
            template_names = []
            for g in glb:
                template_name = g.split('/')[-1].split('.tar')[0]
                if template_name not in pipeline_names:
                    template_names.append(g.split('/')[-1].split('.tar')[0])
            print("pipelines %s ready, %s available."%(', '.join(pipeline_names), ', '.join(template_names)))
            rawpipelineinput = input('which pipeline: ')
            if rawpipelineinput in pipeline_names:
                pipeline = rawpipelineinput
                print("pipeline %s ready for target %s"%(pipeline, target))
                good_input = True
            elif rawpipelineinput in template_names:
                template = rawpipelineinput
                print("extracting template for pipeline %s on target %s"%(template, target))
                commonlib.extract_pipeline(target, template, callerorbenchmark)

    glb = glob.glob('%s/*'%(target_base))
    if '%s/%s'%(target_base, pipeline) in glb:
        pipeline_base = '%s/%s'%(target_base, pipeline)
    else:
        print('pipeline %s not found in target directory %s/'%(pipeline, target))
        sys.exit()

    command_line_args.append(pipeline)

    # build sge submission scripts
    if mode in ['build', 'run']:
        commonlib.build_submission_scripts(target, pipeline)
        maxjobs = 0
        # execute pipeline
        if mode == 'run':
            # check input argument 4: maxjobs
            if not interactive_mode:
                maxjobs = sys.argv[4]
                try:
                    maxjobs = int(sys.argv[4])
                except ValueError:
                    print("maxjobs should be an integer representing the number of threads to spawn. %s provided."%(sys.argv[4]))
                    sys.exit()
            elif interactive_mode:
                good_input = False
                while not good_input:
                    rawmaxjobsinput = input("maximum number of concurrent threads (<200): ")
                    try:
                        maxjobs = rawmaxjobsinput
                        good_input = True
                    except ValueError:
                        print("maxjobs should be an integer representing the number of threads to spawn. %s provided."%(sys.argv[4]))
                        continue
                command_line_args.append(maxjobs)

            print("\nsubmitting jobs for %s"%(' '.join(command_line_args)))
            print("if tictac gets interrupted, just run the following command to restart:")
            print("python tictac.py %s\n"%(' '.join(command_line_args)))
            commonlib.run_submission_scripts(target, pipeline, maxjobs, mode)

        # test mode: attempt to execute the pipeline all steps except actual submission
        elif mode == 'build':
            commonlib.run_submission_scripts(target, pipeline, maxjobs, mode)
            print("check the work/ subdirectory of %s/%s, that the .sh files can be successfully manually submitted."%(target, pipeline))

    if mode in ['check', 'fix']:
        checktype = None
        fixtypes = ['zero_file_size', 'missing_indices']
        fixtype = None
        if mode == 'fix':
            if not interactive_mode:
                fixtype = None
                if len(sys.argv) >= 5:
                    fixtype = sys.argv[4]
                    if fixtype not in fixtypes:
                        print('check type %s not in known check modes %s'%(fixtype, fixtypes))
                else:
                    print("please add a fix type from %s to command line args"%(', '.join(fixtypes)))
                    sys.exit()
            else:
                good_input = False
                while not good_input:
                    rawfixtypeinput = input("fix type (%s): "%(', '.join(fixtypes)))
                    if rawfixtypeinput in fixtypes:
                        fixtype = rawfixtypeinput
                        good_input = True
                    else:
                        print("fix type %s provided. please use types from %s"%(rawfixtypeinput, ', '.join(fixtypes)))
                        continue
                command_line_args.append(fixtype)
        commonlib.check_submission_results(target, pipeline, mode, fixtype)

    if mode == 'train':
        commonlib.train_caller(target, pipeline)

    if mode == 'merge':
        commonlib.merge_vcfs(target, pipeline)

    if mode == 'compare':
        pipelinebase = commonlib._get_pipelinebase(target, pipeline)
        truthvcfsbase = "%s/reference/truthvcfs"%(home_base)
        truthvcftable = {}
        available_truths = glob.glob("%s/*"%(truthvcfsbase))
        a_truths_abbv = []
        for at in available_truths:
            a_truths_abbv.append(at.split('/')[-1])
            truthvcftable[at.split('/')[-1]] = at
        truth_vcf_abbv = None
        if pipeline not in a_truths_abbv:
            # did not find pipeline name in known truths so look for known truths that start the pipeline name
            for ata in a_truths_abbv:
                if pipeline.startswith(ata):
                    truth_vcf_abbv = ata
                    break
            else:
                # did not find any matches so prompt the user for a truth vcf
                good_input = False
                while not good_input:
                    rawtruthinput = input("please specify a truth vcf file location, or add one to reference/truthvcfs: ")
                    if not os.path.isfile(rawtruthinput):
                        print("truth vcf file %s not found. please specify an existing file"%(rawtruthinput))
                        continue
                    else:
                        truth_vcf_abbv = rawtruthinput
                        good_input = True
            # expand the abbreviation to the full location of the truth vcf
            truth_vcf_list = glob.glob('%s/%s/*.vcf.gz'%(truthvcfsbase, truth_vcf_abbv))
            truth_vcf = truth_vcf_list[0]
        else:
            truth_vcf_dir = truthvcftable[pipeline]
            truth_vcf = glob.glob('%s/*.vcf.gz'%(truth_vcf_dir))
            if len(truth_vcf) > 1:
                print("error: %s has more than one '*.vcf.gz' file"%(truth_vcf_dir))
                sys.exit()
            elif len(truth_vcf) == 0:
                print("error: %s has no '*.vcf.gz' file"%(truth_vcf_dir))
                sys.exit()
            truth_vcf = truth_vcf[0]

        print("sending off %s %s %s"%(target, pipeline, truth_vcf))
        commonlib.compare_vcfs_to_truth(target, pipeline, truth_vcf)

    if interactive_mode:
        print("tictac successfully completed job %s"%(' '.join(command_line_args)))
        print("to run the same command from a single command line call, issue the following command at the prompt:")
        print(" $ python tictac.py %s"%(' '.join(command_line_args)))

if __name__ == "__main__":
    main(sys.argv)
