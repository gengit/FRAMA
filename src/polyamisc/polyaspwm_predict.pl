#!/usr/bin/env perl

################################################################################
#
#  polyaspwm_predict.pl
#  score poly(A) signal using PWM prediction
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

our %pas_score = (
  AAAAAA => -1.900,
  AACAAA => -2.600,
  AAGAAA => -2.800,
  AATAAA =>  0.000,
  AATAAC => -3.100,
  AATAAG => -3.000,
  AATAAT => -3.000,
  AATACA => -2.800,
  AATAGA => -3.000,
  AATATA => -2.600,
  AATCAA => -2.800,
  AATGAA => -2.600,
  AATTAA => -1.600,
  ACTAAA => -3.000,
  AGTAAA => -2.900,
  ATTAAA => -1.200,
  ATTTAA => -2.800,
  CATAAA => -3.000,
  GATAAA => -2.600,
  TATAAA => -1.900,
  TTTAAA => -3.100,
  );
our %pas_p;
{ my $i=0;
  foreach (sort { $pas_score{$b}<=>$pas_score{$a} } keys %pas_score) {
    $pas_p{$_} = (++$i)/(4**6);
  }
}

# command line syntax
my $ProgFile = ( split('/',__FILE__) )[-1];
sub usage {
  print "\n";
  print <<END_USAGE;
PURPOSE
 score poly(A) signal using PWM prediction

COMMAND LINE SYNTAX
 $ProgFile [-h]
 $ProgFile [options] seq.fa [seq2.fa [...]] > out.csv

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
    print  "###\n";
    printf "# sequence %s\n", $pSeq->{id};
    my $l=length($pSeq->{seq});
    printf "# seq_length %d\n", $l;
    $pSeq->{seq}=~s/A+$//i;
    $l=length($pSeq->{seq});
    my $bAaa=0;
    for (my $i=0; $i<=$l-6; ++$i) {
      my $s=uc(substr($pSeq->{seq},$i,6));
      if (exists($pas_score{$s})) {
        my $sc=$pas_score{$s};
        my $bOut=1;
        if ($s eq 'AAAAAA') { if ($bAaa) { $bOut=0 } $bAaa=1; } else{ $bAaa=0 }
        if ($bOut) {
          printf "%s\t-\t%d\t%.3f\t%.1e\t%s\n", $pSeq->{id}, $i+1,
            $sc, $pas_p{$s}, $s ;
        }
      }
    }
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
