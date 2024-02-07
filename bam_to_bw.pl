#
use Config::Simple;
use File::Path qw(make_path remove_tree);
##############################################################

($sample) = shift(@ARGV);
load_modules();
$bamDir = "5_nucleosome_free";

$bamFile = "$bamDir/$sample"."_free.bam";
$outBw = "$bamDir/$sample".".bw";
# calc
$cmd = "bamCoverage --bam $bamFile -o $outBw ";
$cmd .= "--binSize 5 --normalizeUsing RPKM ",
   #    "--effectiveGenomeSize 2652783500 ",
#$cmd .=   "--ignoreForNormalization chrX --extendReads";
$cmd .="--ignoreForNormalization chrX";
print "$cmd\n";       
system("$cmd");

sub load_modules{
  do ('/apps/RH7U2/Modules/default/init/perl.pm');
  module('load bamtools');


   return();
}
