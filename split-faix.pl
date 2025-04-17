my $n=$ARGV[1];

open (FILE,$ARGV[0]);
my $totall=0;
my @contigs=();

while(my $line=<FILE>){
  chomp $line;
  my ($c,$l)=split("\t", $line);
  $totall+=$l;
  my $tmp=();
  $tmp->{contig}=$c;
  $tmp->{length}=$l;
  push (@contigs,$tmp);
}

close (FILE);
#print $totall/$n;
my $maxsize=int($totall/$n);
print "MS=$maxsize\n";
my $group=1;
my $acu=0;
open (TT,">pp-$group.targetlist");
if($n==1){
  foreach my $c(@contigs){
      print TT $c->{contig}."\n";
  }
}elsif(scalar(@contigs) == 1){

  foreach my $c(@contigs){
      print TT $c->{contig}."\n";
  }

}else{
foreach my $c(@contigs){
  $acu+=$c->{length};

  if($acu<=$maxsize){
	  #      print join ("\t", $c->{contig}, $c->{length}, $acu, $group)."\n";
      print TT $c->{contig}."\n";
  }
  else{
      $group++;
      $acu=$c->{length};
      #print join ("\t", $c->{contig}, $c->{length}, $acu, $group)."\n";
      close TT;

      open (TT,">pp-$group.targetlist");
      print TT $c->{contig}."\n";
  }


}
}

close TT;
