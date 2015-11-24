#! /usr/bin/perl

use List::Util qw(sum);

$inputfile = $ARGV[0];
$sam_file = $ARGV[1];
$outputfile = $ARGV[2];


open(INPUT, "<$inputfile") || die "Can't open input\n";
open(OUTPUT, ">$outputfile") || die "Can't open output\n";
print OUTPUT "Chromosome\tStart\tEnd\tLog2Ratio\tReads found\tReads expected\tZ-score\n";

$genome_size = 3095693983;
$diploid_size = 5976727269;
$controls_samdir = "/home/peter/Projekte/Plasma_Seq/Kontrollen/Breast/SAM_files";
###############################################################################
sub mean {
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
        print "   Std.-dev: $std\n";
        return $std;
}
###############################################################################
#calculate Z-score of segment
$line = <INPUT>;

while ($line = <INPUT>)
{
  @info = split (/\t/, $line);
  $chr = $info[0];
  $start = $info[1];
  $end = $info[2];
  $length = $end - $start;

  print "Working on Segment: $chr:$start-$end\n";
  open (SAMPLESAM, "<$sam_file") || die "Can't open Samfile\n";
  $count = 0;
  $segment_count = 0;
 
  # get reads on segment 
  while ($samline = <SAMPLESAM>)
  {
    if ($samline =~ m/^@/){
      next;}
    @info = split(/\t/, $samline);
      if ($info[2] eq "*"){
        next;}
      if ($info[4] < 25){
        next;}
      if ($info[1] & 0x4){
        next;}
    $count++;
    if (($info[2] eq $chr) && ($info[3] > $start) && ($info[3] < $end)){
      $segment_count++;}   
  }
  close SAMPLESAM;
  #calculate expected amount of reads by dividing length by genome size and multiply with total aligned reads 
  if (($chr eq "chrX") || ($chr eq "chrY")){
    $expected = ($length / $diploid_size) * $count;}
  else {
    $expected = ($length / $genome_size) * $count;}

  $sample_ratio = $segment_count / $expected;

  print "   Sample processed\treads found: $segment_count\tTotal reads: $count\t ratio: $sample_ratio\n";
  @control_ratios = ();

  #iterate over every control SAM file
  @samfiles = glob("$controls_samdir/*.sam");
  foreach $file (@samfiles)
  {
    open (SAM, "<$file") || die "Can't open Samfile\n";
    $ctrl_count = 0;
    $ctrl_segment_count = 0;
    while ($samline = <SAM>)
    {
      if ($samline =~ m/^@/){
        next;}
      @ctrl_info = split(/\t/, $samline);
      if ($ctrl_info[2] eq "*"){
        next;}
      if ($ctrl_info[4] < 25){
        next;}
      if ($ctrl_info[1] & 0x4){
        next;}

      $ctrl_count++;
      if (($ctrl_info[2] eq $chr) && ($ctrl_info[3] > $start) && ($ctrl_info[3] < $end)){
        $ctrl_segment_count++;}   
    }
    if (($chr eq "chrX") || ($chr eq "chrY")){
      $ctrl_expected = ($length / $diploid_size) * $ctrl_count;}
    else {
      $ctrl_expected = ($length / $genome_size) * $ctrl_count;}
    $ratio = $ctrl_segment_count / $ctrl_expected;

    push(@control_ratios, $ratio);
    close SAM;
    print "   Control $file processed\tReads: $ctrl_segment_count\tExpected: $ctrl_expected\tTotal: $ctrl_count\tRatio: $ratio\n";
  }
  
  $control_mean = &mean(@control_ratios);
  $control_stdev = &stdev(@control_ratios);
  print $control_stdev
  $z_score = ($sample_ratio - $control_mean) / $control_stdev;
  print "   Z-score: $z_score\n";
  $output_line = $line;
  $output_line =~ s/\n//g;
  print OUTPUT "$output_line\t$segment_count\t$expected\t$z_score\n";
  
}

###############################################################################

