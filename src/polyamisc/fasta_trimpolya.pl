#!/usr/bin/env perl
################################################################################
#
#  fasta_trimpolya.pl
#  trim poly(A) stretches in fasta sequence
#
#  copyright (c)
#    FLI Jena, Genome Analysis Group, 2014
#  author
#    Karol Szafranski, szafrans@fli-leibniz.de
#
################################################################################

use strict; #use warnings;  # OK 19700101
use Getopt::Std;
use FileHandle::Unget;

# command line syntax
my $ProgFile = ( split('/',__FILE__) )[-1];
sub usage {
  print "\n";
  print <<END_USAGE;
PURPOSE
 trim poly(A) stretches in fasta sequence

COMMAND LINE SYNTAX
 $ProgFile [-h]
 $ProgFile [options] seq.fa [seq2.fa [...]] > seqout.fa

arguments:
 seq.fa     sequence file in fasta format
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
our (@afseq)=@ARGV;


# loop over input sequences
foreach my $fseq (@afseq) {
  # the open() construct will work with pipes, too (input argument "-")
  open(INSEQ,$fseq) or die "ERROR: cannot read file $fseq";
  my $hSeq=FileHandle::Unget->new(\*INSEQ,'r');
  while (my $pSeq=&read_seq($hSeq)) {
    $pSeq->{seq}=~s/[AN]+$//gi;
    &print_seq($pSeq);
  }
  close(INSEQ);
}

exit;


################################################################################
# functions

# read single sequence entry from fasta file handle
# arg1: file handle for sequence input, object FileHandle::Unget
#       FileHandle object mediates buffered reading from non-physical handle
#       (pipe)
# ret:  hash reference to {id=>$id,descr=>$descr,seq=>$sSeq}
sub read_seq {
  my($hInSeq)=@_;
  while(<$hInSeq>) {
    # initialize sequence entry with fasta header line
    if (m/^>(\S+)(?:\s+(.+))?/) {
      my $pSeq={id=>$1,descr=>$2||''};
      my @seq;
      my $lineo;
      while(defined($_=$lineo=<$hInSeq>) and !m/^>/) {
        chomp;
        push @seq,$_;
      }
      # restore filehandle reading position
      if(length($lineo||'')){ $hInSeq->ungets($lineo) }
      # finish sequence entry
      $pSeq->{seq} = join('',@seq);
      return $pSeq;
    }
  }
  return undef;  # no entry anymore
}

# print single sequence entry in fasta format
# arg1: sequence = hash reference to {id=>$id,descr=>$descr,seq=>$sSeq}
#       referenced data is left unchanged
sub print_seq {
  my($pSeq)=@_;
  printf ">%s %s\n", $pSeq->{id},$pSeq->{descr}||'';
  my $seq=$pSeq->{seq};
  $seq=~s/\s+//g;
  $seq=~s/\S{1,60}/$&\n/g;
  print $seq;
}
