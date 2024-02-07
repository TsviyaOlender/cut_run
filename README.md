# cut_run
pipeline to analyze cut&amp;run NGS data
# This pipeline analyzes cut&run data
The pipeline is written in perl and is designed to work on a cluster, with LSF system.
Dependencies:
fastqc
jdk
R
ngsplot
bowtie2
picard-tools
samtools
bedtools-gnu
python/2.7
bamtools
macs2
###########################################################################
usage:
perl cutrunpipeline.pl samples.txt params_RUNX3.txt
##########################################################################
The script cutrunpipeline is a wrapper that submits the script run_CutRun_V4.pl to the queue. The number of jobs that will be submitted is as the number of 
samples that are specified in the file samples.txt. The file params_RUNX3.txt includes all the run parameters such as the location of the crude reads, genome and more.
To run the script you have to edit line 10, with the full path of the software address in the system,and line 23 with the queue name.
You also have to edit lines 272-282 to refer to the versions of the modules that are installed on your system.
