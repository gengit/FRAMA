#!/usr/local/bin/perl
use strict; #use warnings;  # OK 19700101
use Getopt::Std;

our $colp=0;
our $colsc=1;
our $wsize=4;
our $sc_thresh=0.05;
our $cov_thresh=2400;

# command line usage description
my $ProgFile = ( split('/',__FILE__) )[-1];
sub usage {
  print  "\n";
  print  <<END_USAGE;
PURPOSE
 produce score plot from output of polyaread_sample.pl

COMMAND LINE SYNTAX
 $ProgFile [-h]
 $ProgFile [-v] polyaread.csv

arguments:
 *.csv      report from polyaread_sample.pl

options:
 -c n       table column that holds the score values, default $colsc (computational)
 -d f       consider coverage (/ depth) for normalization of polyA read frequency.
            Give file with data produces by script polyacov_coverage.pl
 -h         print command line usage description and exit
END_USAGE
  print  "\n";
  exit 1;
}

# command-line interface
our %mopt;
getopts('c:d:h',\%mopt) or &usage();
if ($mopt{h} or int(@ARGV)<1) { &usage() }
if ($mopt{c}) { $colsc=$mopt{c} }
our ($fin)=@ARGV;


our %diccov;
if ($mopt{d}) {
  my $s=undef;
  open(IN,$mopt{d}) or die "ERROR: cannot read file $mopt{d}";
  while (<IN>) {
    if (m/^#/){ next }
    chomp;
    my @v=split /\t/;
    $diccov{$v[0]}[$v[1]]=$v[2];
  }
  close(IN);
}

open(IN,$fin) or die "ERROR: cannot read file $fin";
our $s=undef;
our $p-1;
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
        if (m/^#\s*sequence\s+(\S+)/){ $s=$1 }
      }
      else{ last }
    }
    redo;
  }
  chomp;
  # read occurrence of polyA reads
  my @v=split/\t/;
  if ($mopt{d}) {
    my $mcov=&mean ( @{$diccov{$s}}[&max($v[$colp]-10,1)..&max($v[$colp],10)] );
    if ($mcov>$cov_thresh) {
      $v[$colsc] /= $mcov/$cov_thresh;
    }
  }
  # next item in current sequence
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

exit;


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
       * $pCand->[$main::colsc];
  }
  if ($sc>$sc_thresh){
    printf "%s\t%d\t%.3f\t%s%s\n", $main::s, $p, $sc,
      join(',',map{ $_->[$main::colp] } @main::cand), $main::mopt{d}?sprintf("\t%d",$main::diccov{$main::s}[$p]):'';
  }
}
