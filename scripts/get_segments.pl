#! /usr/bin/perl

$inputfile = $ARGV[0];
$outputfile = $ARGV[1];


open(INPUT, "<$inputfile") || die "Can't open input\n";
open(OUTPUT, ">$outputfile") || die "Can't open output\n";
print OUTPUT "Chromosome\tStart\tEnd\tLog2Ratio\n";

$lastseg = 0;
$lastpos = 0;
$firstpos = 0;
$last_chr = "";

$line = <INPUT>;

while ($line = <INPUT>)
{
  @info = split (/\t/, $line);
  $seg = $info[6];
  $seg =~ s/\n//g;
  #print "$lastseg $seg\t$info[6]";
  if ($seg eq $lastseg){
    $lastpos = $info[3];
    $lastseg = $seg;
    next;}
  elsif (($lastseg != 0) && ($seg != $lastseg))
  {
    if ($last_chr eq ("chr".$info[1])){
      $lastpos = $info[3]-1;}
    if ($last_chr eq "chr23"){
      $last_chr = "chrX";}
    elsif($last_chr eq "chr24"){
      $last_chr = "chrY";}
    print OUTPUT "$last_chr\t$firstpos\t$lastpos\t$lastseg\n";
    $firstpos = $info[3];
    $last_chr = "chr".$info[1];
    $lastseg = $seg;
  }
  elsif (($lastseg == 0) && ($seg ne $lastseg))
  {
    $firstpos = $info[3];
    $lastseg = $seg; 
    $last_chr = "chr".$info[1];
  }    
}

if ($last_chr eq "chr23"){
  $last_chr = "chrX";}
elsif($last_chr eq "chr24"){
  $last_chr = "chrY";}
print OUTPUT "$last_chr\t$firstpos\t$lastpos\t$lastseg\n";
close OUTPUT;
close INPUT;


