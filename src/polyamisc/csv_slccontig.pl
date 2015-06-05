#! /usr/local/bin/perl
################################################################################
#
#  csv_slccontig.pl
#  select contig passage from polya prediction csv stream
#
#  copyright (c)
#    FLI Jena, Genome Analysis Group, 2014
#  author
#    Karol Szafranski, szafrans@fli-leibniz.de
#
################################################################################

use strict; #use warnings;  # OK 19700101
use Getopt::Std;

# command line syntax
my $ProgFile = ( split('/',__FILE__) )[-1];
sub usage {
  print "\n";
  print <<END_USAGE;
PURPOSE
 select contig passage (1st occurrence) from polya prediction csv stream

COMMAND LINE SYNTAX
 $ProgFile [-h]
 $ProgFile [options] contigid polya*.csv [polya*2.csv [...]] > out.csv

arguments:
 contigid   contig identifier for selection
 polya*.csv polya prediction data in *.csv format
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
if ($mopt{h} or int(@ARGV)<2) { &usage() }
our ($ctg,@afin)=@ARGV;


# loop over input sequences
foreach my $fin (@afin) {
  # the open() construct will work with pipes, too (input argument "-")
  open(IN,$fin) or die "ERROR: cannot read file $fin";
  while (<IN>) {
    if (m/^###/) {
      $_=<IN>;
      if (m/^#\s*(?:sequence\s+)?(\S+)/ and $1 eq $ctg) {
        print "###\n";
        print;
        $_=<IN>; print;
        while (<IN>) {
          if (m/^#/) { last }
          else { print }
        }
      }
    }
  }
  close(IN);
}
