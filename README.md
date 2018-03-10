Tictac automates variant calling. It segment a bam into intervals for parallelization, runs callers, concatenates output vcfs, and compares to truth. It requires SGE, common bioinformatics tools like samtools, bcftools, bedtools, and picard, and some variant callers.

I have templates for pipelines that do variant calling, bam-segmenting, and pileups, all of which I will add soon. 

Install by cloning the repo. Add environmental variables for TICTAC_HOME (tictac's head directory) and TICTAC_QUEUE. Tictac relies heavily on SGE. That should be all you need to get started.
