#!
my($input) = @ARGV;
open(IN,"$input");
do{
  $line = <IN>;
}until($line =~/^chr/);
while($line = <IN>){
  chomp($line);
  (@data) = split(/\t/,$line);
  print "$data[0]\t$data[1]\t$data[2]\t$data[9]\n";
}