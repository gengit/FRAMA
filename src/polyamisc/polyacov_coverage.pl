#!/usr/bin/env perl
################################################################################
#
#  polyacov_coverage.pl
#  score poly(A) signal using coverage information
#
#  copyright (c)
#    FLI Jena, Genome Analysis Group, 2014
#  author
#    Karol Szafranski, szafrans@fli-leibniz.de
#
################################################################################

use strict;
use Getopt::Std;
use FileHandle::Unget;

our $log2=log(2);

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
    $seqid=$v[0];
    print  "###\n";
    printf "# sequence %s\n", $seqid;
  }
  printf "%s\t%d\t%.3f\n", $v[0], $v[1], log($v[2]+0.99)/$log2;
}
close(IN);
