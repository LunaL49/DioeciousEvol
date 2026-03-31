This folder contains all the files associated with running the bioinformatic pipeline.

To run the bioinformatics pipeline:<br>
Step 1: Run preprocess.sh on paired-end raw Illumina sequencing files. This outputs 20 reverse demultiplexed files which need to be processed in pairs in step 2.<br>
Step 2: Run either ben1_call_variants.sh or pod2_call_variants.sh with the name of the pair as command line argument, e.g. "R2" or "R6". <br>

See bash scripts for more information on calling syntax and input/output files. 
