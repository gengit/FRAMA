#! /usr/local/bin/perl
################################################################################
#
#  polya_coverage.pl
#  score poly(A) signal using coverage information
#
#  copyright (c)
#    FLI Jena, Genome Analysis Group, 2014
#  author
#    Karol Szafranski, szafrans@fli-leibniz.de
#
################################################################################

use strict; #use warnings;  # OK 20140630
use Getopt::Std;
use FileHandle::Unget;

our $wsize=5;
our $sc_thresh=0.05;
our $colp=0;
our $colsc=1;
our $pk_shift=-3;  # ad hoc 0

# command line syntax
my $ProgFile = ( split('/',__FILE__) )[-1];
sub usage {
  print "\n";
  print <<END_USAGE;
PURPOSE
 score poly(A) signal using PWM prediction

COMMAND LINE SYNTAX
 $ProgFile [-h]
 $ProgFile [options] bam2coverage.csv > out.csv

arguments:
 bam2coverage.csv
            coverage data in *.csv format, as produced by bamtools -f coverage
            input argument "-" evaluates to STDIN

options:
 -h         show this command-line syntax description and exit
END_USAGE
  print "\n";
  exit 1;
}

# command-line interface
our %mopt;
getopts('hx:',\%mopt) or &usage();
if ($mopt{h} or int(@ARGV)<1) { &usage() }
our ($fin)=@ARGV;


open(IN,$fin) or die "ERROR: cannot read file $fin";
our $s=undef;
our $p = -1;
our @cand;
while (<IN>) {
  if (m/^#/){
    # close previous sequence
    if (length($s)){
      foreach ( ; int(@cand); ++$p){ &score_eval($p) }
      $p=-1;
      @cand=();
    }
    print;
    undef $s;
    # initialize for next sequence
    while (<IN>) {
      if (m/^#/){
        print;
        if (m/^#\s*(\S+)/){ $s=$1 }
      }
      else{ last }
    }
    redo;
  }
  chomp;
  my @v=split/\t/;
  # next item in current sequence
  $v[$colsc]=&score_func($v[$colsc]/100);
  if (int(@cand)) {
    my $pn=&min($v[$colp]-2*$wsize,$cand[-1][$colp]+2*$wsize);
    push @cand, \@v;
    foreach ( ; $p<$pn; ++$p){ &score_eval($p) }
  }
  else{
    $p=&max(($p>0)?$p:(),$v[$colp]-2*$wsize);
    my $pn=$v[$colp]-2*$wsize;
    push @cand, \@v;
    foreach ( ; $p<$pn; ++$p){ &score_eval($p) }
  }
}

foreach ( ; int(@cand); ++$p){ &score_eval($p) }


################################################################################
# misc. functions

# minimum of array of values
sub min {
  my ($min,@v) = grep{ defined($_) } @_;
  foreach (@v) {
    $min = ($min<$_) ? $min : $_;
  }
  return $min;
}

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


# calculate mean value from array of values
sub score_func {
  my ($v)=@_;
  my @dicfunc=(
    [0.000,0.000],
    [0.590,0.0005],
    [0.600,0.001],
    [0.610,0.002],
    [0.620,0.003],
    [0.630,0.004],
    [0.640,0.005],
    [0.650,0.006],
    [0.660,0.009],
    [0.670,0.013],
    [0.680,0.021],
    [0.690,0.034],
    [0.700,0.058],
    [0.710,0.103],
    [0.720,0.184],
    [0.730,0.320],
    [0.740,0.510],
    [0.750,0.708],
    [0.760,0.856],
    [0.770,0.939],
    [0.780,0.976],
    [0.790,0.991],
    [0.800,0.996],
    [0.810,0.997],
    [0.820,0.998],
    [0.830,0.998],
    [0.840,0.998],
    [0.850,0.998],
    [0.860,0.998],
    [0.870,0.998],
    [0.880,0.998],
    [0.890,0.998],
    [0.900,0.999],
    [0.910,0.999],
    [0.920,0.999],
    [0.930,0.999],
    [0.940,0.999],
    [0.950,0.999],
    [0.960,0.999],
    [0.970,0.999],
    [0.980,0.999],
    [0.990,0.999],
    [1.000,1.000],
    [1.001,1.000], 
    );
  # locate x low-bound in dictionary
  my $i; for ($i=$#dicfunc-1; $v<$dicfunc[$i][0]; --$i) { }
  # interpolate
  my $ret = $dicfunc[$i][1] + ($dicfunc[$i+1][1]-$dicfunc[$i][1]) * ($v-$dicfunc[$i][0])/($dicfunc[$i+1][0]-$dicfunc[$i][0]);
  #print "$dicfunc[$i][1] + ($dicfunc[$i+1][1]-$dicfunc[$i][1]) * ($v-$dicfunc[$i][0])/($dicfunc[$i+1][0]-$dicfunc[$i][0]) = $ret\n";
  return $ret;
}

sub gaussval {
  my ($my,$sig,$x) = @_;
  if ($sig<=0) { return undef }
  # calculate value
  my $g = 0.3989422804 * exp -(($x-$my)**2 / (2 * $sig**2));
  return $g;
}


sub score_eval {
  my ($p) = @_;
  while (int(@main::cand) and $main::cand[0][$main::colp]+2*$main::wsize<$p) { shift @main::cand }
  my $sc=0;
  foreach my $pCand (@main::cand) {
    $sc += &gaussval($pCand->[$main::colp],$main::wsize,$p)
       * 2**$pCand->[$main::colsc];
  }
  if ($sc>$sc_thresh){
    printf "%s\t%d\t%.3f\t%s\n", $main::s, $p+$main::pk_shift, $sc,
      join(',',map{ $_->[$main::colp] } @main::cand);
  }
}
