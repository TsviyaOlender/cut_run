#!
##############################
use Config::Simple;
use File::Path qw(make_path remove_tree);
##############################################################
# the script submits samples to the program run_ATAC_Ts_V1.pl
# the input is a file with sample names
#################################################################
my($samplesFile,$paramFile,$queue,$mem) = (@ARGV);
$progD = "/home/labs/home_dir/Cut_Run";
open(IN,"$samplesFile")|| die "can not locate file with sample names on disk\n";
@samples = <IN>;
chomp(@samples);
close(IN);

#set up default memory
if($mem > 0){
}else{
  $mem=8000;
}
if($queue =~/[a-z]/){
}else{
  $queue = 'molgen-short';
}
## prepare run folders
make_folders();

foreach $sample (@samples){
  $logF = "$sample"."_runlog.txt";
  $errF = "$sample"."_err.txt";
  
  $cmd ="bsub -q $queue -n 8 -o $logF -e $errF -R \"span[hosts=1]\" -R \"rusage[mem=$mem]\" perl $progD/run_CutRun_V4.pl $sample $paramFile $progD";
  print "$cmd\n";
  system("$cmd");
}
####################################################################################
sub make_folders{
  make_path("1_fastq");
  make_path("2_fastqc");
  make_path("3_align");
  make_path("4_plots");
  make_path("5_nucleosome_free");
  make_path("6_MACS_2");
  make_path("7_TSS");
 return();
}