#! /usr/bin/perl

use List::Util qw(sum);

$segmentsfile=$ARGV[0];
$bincount_file=$ARGV[1];
$outfile = $ARGV[2];
$control_dir=$ARGV[3];#"/home/peter/Projekte/Plasma_Seq/Kontrollen/Breast/GC_corrected_bincounts";

open (OUTFILE, ">$outfile") || die "Can't open Output file\n";
print OUTFILE "Chromosome\tStart\tEnd\tRatio-Sum\tZ-score\n";
#######################################################################################
sub mean
{ 
  return sum(@_)/@_;
}
###############################################################################
sub stdev{
        @data = @_;
        if(scalar(@data) == 1){
                print "only one data point\n";
                return 0;
        }
        $average = &mean(@data);
        my $sqtotal = 0;
        foreach $datapoint (@data) 
        {
          #print "Datapoint $datapoint\n";
          $sqtotal += ($average-$datapoint) ** 2;
        }
        my $std = ($sqtotal / (scalar(@data))) ** 0.5;
        #print "   Std.-dev: $std\n";
        return $std;
}
#######################################################################################
sub getSumRatio
{
  my $filename  = $_[0];
  my $sub_chr   = $_[1]; 
  my $sub_start = $_[2]; 
  my $sub_end   = $_[3]; 
  my $return_sum = 0;
  open (FILE, "<$filename") || die "Can't open $filename\n";
  $sub_line = <FILE>;
  while ($sub_line = <FILE>)
  {
    @file_info = split(/\t/, $sub_line);
    if ($file_info[0] eq "23"){
      $actual_chr = "chrX";}
    elsif ($file_info[0] eq "24"){
      $actual_chr = "chrY";}
    else {
      $actual_chr = "chr".$file_info[0];}
 
    if (($sub_chr eq $actual_chr) && ($sub_start <= $file_info[1]) && ($sub_end >= $file_info[1]))
    {
      $sum =  $file_info[6];
      $sum =~ s/\n//g;
      $return_sum += $sum;
    }
  }
  close FILE;
  return $return_sum;
  
}
#######################################################################################
@control_files=glob("$control_dir/*.bincounts");

open(SEGMENTS, "<$segmentsfile") || die "Can't open segments\n";
$line=<SEGMENTS>;
while ($line=<SEGMENTS>)
{
  @info = split("\t", $line);
  $chromosome = $info[0];
  $start = $info[1];
  $end = $info[2];
 
  $sumratio = getSumRatio($bincount_file, $chromosome, $start, $end);


  my @controls_sum;  
  foreach $control_file (@control_files)
  {
    push(@controls_sum, getSumRatio($control_file, $chromosome, $start, $end));
  }

  $average = mean(@controls_sum);
  $stdev = stdev(@controls_sum);
  $z_score_Segment = ($sumratio - $average) / $stdev;
  print OUTFILE "$chromosome\t$start\t$end\t$sumratio\t$z_score_Segment\n";
}

close OUTFILE;
