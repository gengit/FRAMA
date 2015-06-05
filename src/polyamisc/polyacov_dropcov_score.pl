#! /usr/local/bin/perl
################################################################################
#
#  polyacov_dropcov_score.pl
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

our $log2=log(2);
our $cov1_flank=30;
our $cov1_blank=5;
our $cov2_flank=10;
our $cov2_blank=2;
our $pseudo=1.5;

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


# input csv table; produced by bamtools coverage
# the open() construct will work with pipes, too (input argument "-")
open(IN,$fin) or die "ERROR: cannot read file $fin";
our $pBufferW;
our $seqid;
while (<IN>) {
  chomp;
  my @v=split/\t/;
  if ($v[0] ne $seqid) {
    $pBufferW=CovBufferWide->new();
    $seqid=$v[0];
    print  "###\n";
    printf "# sequence %s\n", $seqid;
    for (my $i=1.6*$main::cov1_flank+$main::cov1_blank; $i>0; --$i) {
      $pBufferW->push_pos($v[1]-$i,0);
    }
  }
  $pBufferW->push_pos($v[1],$v[2]);
}
close(IN);


################################################################################
# misc. functions

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


################################################################################
package CovBufferWide;

sub new {
  my (undef,$f) = @_;
  my $this = {buffer=>[],buffcorr=>CovBufferSharp->new()};
  bless $this;
}

sub push_pos {
  my ($this,$p,$cov) = @_;
  push @{$this->{buffer}}, [$p,$cov+1];
  if (int(@{$this->{buffer}}) > 2*$main::cov1_flank+$main::cov1_blank) {
    shift @{$this->{buffer}};
    my $covl=&main::mean(map{ $_->[1]+$main::pseudo } @{$this->{buffer}}[0..$main::cov1_flank]);
    my $covr=&main::mean(map{ $_->[1]+$main::pseudo } @{$this->{buffer}}[-$main::cov1_flank..-1]);
    my $po=$this->{buffer}[$#{$this->{buffer}}*0.5][0];
    if ($po>0) {
      my $corr=$this->{buffcorr}->query_pos($po);
      printf "%s\t%d\t%.3f\n", $main::seqid, $po, log($covl/$covr)/$main::log2-$corr*2;
    }
  }
}

sub DESTROY {
  my ($this) = @_;
  for (my $i=0; $i<(1.6*$main::cov1_flank+$main::cov1_blank); ++$i) {
    $this->push_pos($this->{buffer}[-1][0]+1,0);
  }
  return;
}


################################################################################
package CovBufferSharp;

sub new {
  my (undef,$f) = @_;
  my $this = {buffer=>[],history=>[]};
  bless $this;
}

sub push_pos {
  my ($this,$p,$cov) = @_;
  push @{$this->{buffer}}, [$p,$cov+1];
  if (int(@{$this->{buffer}}) > 2*$main::cov2_flank+$main::cov2_blank) {
    shift @{$this->{buffer}};
    my $covl=&main::mean(map{ $_->[1]+$main::pseudo } @{$this->{buffer}}[0..$main::cov2_flank]);
    my $covr=&main::mean(map{ $_->[1]+$main::pseudo } @{$this->{buffer}}[-$main::cov2_flank..-1]);
    my $po=$this->{buffer}[$#{$this->{buffer}}*0.5][0];
    if ($po+3>0) {  # index shifted +3, to keep buffer for value smoothening
      $this->{history}[$po+3] = log($covl/$covr)/$main::log2;
    }
  }
}

sub query_pos {
  my ($this,$p) = @_;
  $p+=3;  # index shifted +3
  # smoothening by triangle profile, size +3
  my $v=0; my $f=0; my $sz=3;
  for (my $i=0; $i<=$sz; ++$i) {
    my $lf = ($sz+1-$i) / 10;
    $f += $lf;
    $v += $lf*$this->{history}[$p-$i];
    if ($i>0) {
      $f += $lf;
      $v += $lf*$this->{history}[$p+$i];
    }
  }
  return $v/$f;
}

sub DESTROY {
  my ($this) = @_;
  for (my $i=0; $i<(1.6*$main::cov2_flank+$main::cov2_blank); ++$i) {
    $this->push_pos($this->{buffer}[-1][0]+1,0);
  }
  return;
}
