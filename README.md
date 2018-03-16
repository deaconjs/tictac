Tictac automates variant calling. It segment a bam into intervals to parallelize callers. It then runs callers and checks, fixes, and concatenates output vcfs. It requires SGE, common bioinformatics tools like samtools, bcftools, bedtools, and picard, and some variant callers. Few to lots of cores is great.

Templates for variant callers will be available through the Thousand Variant Callers [Repo](https://github.com/deaconjs/ThousandVariantCallersRepo). Templates for utility pipelines that segment bams into intervals, and that run samtools pileup on those intervals will be available here. 

Clone the repo to install. Add environmental variables for TICTAC_HOME (tictac directory) and TICTAC_QUEUE (SGE queue). 

The tictac.py file is run alone to enter interactive mode. There a set of commands are elicited from the user to select a target genome and pipeline. It helps write the pipeline executable scipt with templates. The same commands 

