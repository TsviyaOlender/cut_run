#!
use Bio::Seq;
use Bio::SeqIO;

my($sample,$genome) = @ARGV;
$bedfile=$sample.".bed";
$tmpF = $sample."TMP.txt";
$outFile = $sample.".RBS_TFBS.bed";
$outFile_RBS = $sample.".TFBS_runx3_f.txt";
#$genome = "/home/labs/olenderlab/lvzvia/WIS_exp/genome/mm10.fa";
#$genome = "/shareDB/genomes/Mus_musculus/UCSC/mm9/mm9.fa";
#$genome = "/shareDB/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa";
#$genome = "/shareDB/iGenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa";
load_peaks();
#load genome
$in  = Bio::SeqIO->new(-file => "$genome" ,
                               -format => 'Fasta');
@RUNXsites=();
@RUNXcSites=();
open(OUT,">$tmpF");
while($chrSeq= $in->next_seq() ) {
  $chr=$chrSeq->id;
  foreach my $peak (keys %{$peaks{$chr}}){
      $candidate_seq = $chrSeq->subseq($peaks{$chr}{$peak}{start},$peaks{$chr}{$peak}{end});
      $candidate_seq_rev = get_revcom($candidate_seq);
      
      #UC
      $candidate_seq=uc($candidate_seq);
      $candidate_seq_rev =uc($candidate_seq_rev);
     #search Runx3
     $Pstrand='plus';
     @RUNXsites=search_for_motif1($candidate_seq,$peaks{$chr}{$peak}{start},$Pstrand,$peak);
      $Pstrand='minus';
     #@RUNXcSites=search_for_motif1($candidate_seq_rev,$peaks{$chr}{$peak}{end},$Pstrand); 
 } 
}
close(OUT);

# analyze
open(IN,"$tmpF");
open(OUT1,">$outFile");
%peak_count = ();
while($line = <IN>){
  chomp($line);
  ($chr,$start,$end,$name,$score,$strand,$peakName)=split(/\t/,$line);
  print OUT1 "$chr\t$start\t$end\t$name\t$score\t$strand\n";
  $peak_count{$peakName}++;
  
}
close(IN);
close(OUT1);

open(OUT2,">$outFile_RBS");
foreach $peak (keys %peak_count){
  if($peak_count{$peak} > 0){
    print OUT2 "$peak\t$peak_count{$peak}\n";
  }else{
    print OUT2 "$peak\t0\n";
  }
}
close(OUT2);

system("rm $tmpF");
################################################################
sub search_for_motif1{
  my($peak_seq,$peakLocation,$flag,$peakName) = @_;
#5’-TG(T/C)GG(T/C)=> 
#(A/G)CC(A/G)CA
 # @motifs = ('ATAAATAAT');
#  @motifs =('ATAAATAAT','ATAAATTAT','ATAATTAAT','ATAATTTAT','ATTAATAAT','ATTAATTAT','ATTATTAAT','ATTATTTAT');

  @motifs = ('TGTGGT','TGCGGT','TGTGGC','TGCGGC','ACCACA','ACCGCA','GCCACA','GCCGCA');
  $seqL=length($peak_seq);
  $motifL = 6;  
  
  for($i=0;$i<$seqL-$motifL+1;$i=$i+1){
    foreach $motif (@motifs){
      my($MySeq)=substr($peak_seq,$i,$motifL);
      
      if($MySeq eq $motif){
      #if($MySeq =~/([AGC][AGC][CG]TCAT[AT]AATTAT)/){
    
        if($flag eq 'plus'){
          $motifLoc = $peakLocation+$i;
          $motifEnd = $motifLoc+$motifL-1;
          $peakStr="+";
        }elsif($flag eq 'minus'){
          $motifLoc = $peakLocation-$i-($motifL-1);
          $motifEnd = $motifLoc+$motifL-1;
          $peakStr="-";
        }
        print OUT "$chr\t$motifLoc\t$motifEnd\t$MySeq\t1\t$peakStr\t$peakName\n";
       
      }
    } 
  }
  $i=0;

  return();
}

sub get_revcom{
  my($inSeq) = @_;
  my($out_revcom) = '';
  $seqobj = Bio::Seq->new( -seq => $inSeq);
  $revcomobj = $seqobj->revcom();
  $out_revcom=$revcomobj->seq();

  return($out_revcom);
}


sub load_peaks{
  open(IN,"$bedfile")|| warn "can not loacte bedfile\n";
  while($line = <IN>){
    chomp($line);
    (@data) = split(/\t/,$line);
    $peakName = $data[3];
    $peaks{$data[0]}{$peakName}{start} = $data[1];
    $peaks{$data[0]}{$peakName}{end} = $data[2];
  
  }
  close(IN);
  return();
}