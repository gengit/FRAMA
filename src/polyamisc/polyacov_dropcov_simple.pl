#! /usr/local/bin/perl
################################################################################
#
#  polyacov_dropcov_simple.pl
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
our $cov_flank=30;
our $cov_blank=2;
our $pseudo=1.5;
our $pk_shift=4;  # ad hoc 0

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
our $pBuffer;
our $seqid;
while (<IN>) {
  chomp;
  my @v=split/\t/;
  if ($v[0] ne $seqid) {
    $pBuffer=CovBuffer->new();
    $seqid=$v[0];
    print  "###\n";
    printf "# sequence %s\n", $seqid;
    for (my $i=1.5*$main::cov_flank+$main::cov_blank; $i>0; --$i) {
      $pBuffer->push_pos($v[1]-$i,0);
    }
  }
  $pBuffer->push_pos($v[1],$v[2]);
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
package CovBuffer;

sub new {
  my (undef,$f) = @_;
  my $this = [];
  bless $this;
}

sub push_pos {
  my ($this,$p,$cov) = @_;
  # filter extreme jumps in coverage (2**0.9 or 2**-0.9)
  if (int(@{$this}) and abs(log($this->[-1][1]/($cov+1))/$main::log2)>0.9) {
    my $f = ($this->[-1][1]>$cov)? 2**(0.9-log($this->[-1][1]/($cov+1))/$main::log2) : 2**(-0.9-log($this->[-1][1]/($cov+1))/$main::log2);
    for (my $i=0; $i<int(@{$this}); ++$i) {
      $this->[$i][1]*=$f;
    }
  }
  # enter next value to buffer
  push @{$this}, [$p,$cov+1];
  # calculate window center on filled buffer
  if (int(@{$this}) > 2*$main::cov_flank+$main::cov_blank) {
    shift @{$this};
    my $covl=&main::mean(map{ $_->[1]+$main::pseudo } @{$this}[0..$main::cov_flank]);
    my $covr=&main::mean(map{ $_->[1]+$main::pseudo } @{$this}[-$main::cov_flank..-1]);
    printf "%s\t%d\t%.3f\n", $main::seqid,
      $this->[$#{$this}*0.5][0]+$main::pk_shift, log($covl/$covr)/$main::log2;
  }
}

sub DESTROY {
  my ($this) = @_;
  for (my $i=0; $i<1.5*$main::cov_flank; ++$i) {
    $this->push_pos($this->[-1][0]+1,0);
  }
  return;
}
