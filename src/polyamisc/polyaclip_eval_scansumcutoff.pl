#!/usr/bin/env perl
use strict;
use Getopt::Std;

our $scheme=0;

# command line usage description
my $ProgFile = ( split('/',__FILE__) )[-1];
sub usage {
  print  "\n";
  print  <<END_USAGE;
PURPOSE

COMMAND LINE SYNTAX
 $ProgFile [-h]
 $ProgFile [options] polya_genes.csv polya_summary.csv [polya_summary2.csv [...]]

arguments:
 polya_genes.csv
            gene info dictionary produced by polyaclip_eval_geneprop.pl
 polya_summary.csv
            sum scoring produced by polyaclip_sumscore.pl

options:
 -1         extra reward number of clippings
 -h         print command line usage description and exit
END_USAGE
  print  "\n";
  exit 1;
}

# command-line interface
our %mopt;
getopts('1hx:',\%mopt) or &usage();
if ($mopt{h} or int(@ARGV)<2) { &usage() }
if ($mopt{1}) { $scheme=1 }


# load dictionary of optimality
my $fdic=shift();
our %dic;
open(IN,$fdic) or die "ERROR: cannot read file $fdic\n";
while (<IN>) {
  if (m/^#/){ next }
  chomp;
  my @v=split/\t/;
  $dic{$v[0]} = [@v];
}
close(IN);

# determine winner positions from sumscore.csv (per gene)
our %pred;
foreach my $fin (@ARGV) {
  my ($g) = ($fin=~m/.*\/(\w\S+?)_/);
  my @cand = (-1,0);
  open(IN,$fin) or die "ERROR: cannot read file $fin\n";
  while (<IN>) {
    if (m/^#/){ next }
    my @v=split/\t/;
    if ($v[2]>$cand[1]) { @cand=($v[1],$v[2]) }
  }
  close(IN);
  $pred{$g} = [@cand];
}

print  "#thresh\tpenalty[k]\tnclipok\tnclip\tnclipshort\tnclipbad(clipbad)\tXXX\n";
my @aThr;
for (my $thr=0.64; $thr<1.50; $thr+=0.02) { push @aThr,$thr }
for (my $thr=0.96; $thr<1.14; $thr+=0.005) { push @aThr,$thr }
my $thrprev;
foreach my $thr (sort{ $a<=>$b } @aThr) {
  if ($thrprev and abs($thr-$thrprev)<0.0001) { next }
  my $penal=0;
  my $nc=0; my $csh=0; my $cok=0; my $cbad=0; my @cval; my @gbad;
  foreach my $g (sort keys %dic) {

    # winner position, accounting for score threshold (default: leave 3end)
    my $p = ($pred{$g}[1]>$thr)? $pred{$g}[0] : $dic{$g}[1]+1;

    # distance-based penalty for clip prediction
    $penal += ($p<$dic{$g}[4])? 1.5*log($dic{$g}[4]-$p) : log(&max($p-$dic{$g}[4],0.5));

    # classification of clip prediction distances
    if ($pred{$g}[1]>$thr) { ++$nc }
    if ($pred{$g}[1]>$thr and $p<$dic{$g}[1]-10) {
      ++$csh;
      if ($p<$dic{$g}[4]-40) { ++$cbad; push @gbad,$g; }
      push @cval, $dic{$g}[1]-$p;
    }
    if (abs($p-$dic{$g}[4])<=25) { ++$cok }
  }

  # total score, report
  my $total = $scheme? $cok * ($nc**1.8/40) / $penal : $cok * ($nc**1.1/1.5) / $penal;
  printf "%.3f\t%.2f\t%d\t%d\t%d\t%d(%s)\t%.1f\n", $thr, $total,
    $cok,$nc,$csh,$cbad,($cbad<15)?join(',',@gbad):'multiple', int(@cval)?&mean(@cval):0;
  $thrprev=$thr;
}


################################################################################
# misc. functions

# maximum of array of values
sub max {
  my ($max,@v) = grep{ defined($_) } @_;
  foreach (@v) {
    $max = ($max>$_) ? $max : $_;
  }
  return $max;
}

# sum of values
sub sum {
  my ($sum,@v) = grep{ defined($_) } @_;
  foreach (@v) { $sum+=$_ }
  return $sum;
}

# calculate mean value from array of values
sub mean {
  my (@data)=@_;
  my $n=int(@data);
  unless ($n) { die sprintf "mean. ERROR: missing the data!\n" }
  return &sum(@data) / $n;
}
