

This is the home of tictac.

Tictac automates fast variant calling. Pipelines are ok for target, exome, whole genome, somatic, germline, pretty much any variant calling. Tictac runs on local machines for now, but soon will be cloud-optional.

Tictac speeds up processing by segmenting a sequencing file into regions, which can easily be run in parallel. 20 cores means 20x faster calls. Tictac runs the callers then checks and fixes the output. It then concatenates output regional calls into one vcf. It requires SGE, common bioinformatics tools like samtools, bcftools, bedtools, and picard, and some variant callers. 

Templates for variant callers will be available through the Thousand Variant Callers [Repo](https://github.com/deaconjs/ThousandVariantCallersRepo). Templates for utility pipelines that segment bams into intervals, and that run samtools pileup on those intervals will be available here. 

Clone the repo to install, then add environmental variables for TICTAC_HOME (tictac directory) and TICTAC_QUEUE (SGE queue). 

The tictac.py file is run alone to enter interactive mode. There a set of commands are elicited from the user to select a target genome and pipeline. It helps write the pipeline executable scipt with templates. The same commands 

