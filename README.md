# cut_run
pipeline to analyze cut&amp;run NGS data<br>
# This pipeline analyzes cut&run data<br>
The pipeline is written in perl and is designed to work on a cluster, with LSF system.<br>
Dependencies:<br>
fastqc<br>
jdk<br>
R<br>
ngsplot<br>
bowtie2<br>
picard-tools<br>
samtools<br>
bedtools-gnu<br>
python/2.7<br>
bamtools<br>
macs2<br>
###########################################################################<br>
usage:<br>
perl cutrunpipeline.pl samples.txt params_RUNX3.txt<br>
##########################################################################<br>
The script cutrunpipeline is a wrapper that submits the script run_CutRun_V4.pl to the queue. The number of jobs that will be submitted is as the number of <br>
samples that are specified in the file samples.txt. The file params_RUNX3.txt includes all the run parameters such as the location of the crude reads, genome and more.<br>
To run the script you have to edit line 10, with the full path of the software address in the system,and line 23 with the queue name.<br>
You also have to edit lines 272-282 to refer to the versions of the modules that are installed on your system.<br>
