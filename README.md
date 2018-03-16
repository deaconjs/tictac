Tictac automates variant calling. It segment a bam into intervals for parallelization, runs callers, concatenates output vcfs, and compares to truth. It requires SGE, common bioinformatics tools like samtools, bcftools, bedtools, and picard, and some variant callers.

Templates for bam-segmenting, and mpileups will be available here. Templates for variant callers will be available through the Thousand Variant Callers [Repo](https://github.com/deaconjs/ThousandVariantCallersRepo).

Clone the repo to install. Add environmental variables for TICTAC_HOME (tictac directory) and TICTAC_QUEUE (SGE queue). 

The tictac.py file can be run alone to enter interactive mode, wherein a set of commands selects a target genome and helps write the executable scipt. Then the same set expressed as command line arguments re-runs that pipeline. Available commands are listed at input.

