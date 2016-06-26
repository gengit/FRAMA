#!/usr/local/bin/perl
################################################################################
#
#  FRAMA package
#  find_readlist.pl
#  in Trinity data directory, find read list files associated with final contigs
#
#  copyright (c)
#    FLI Jena, Genome Analysis Group, 2016
#  author
#    Karol Szafranski, szafrans@fli-leibniz.de
#
################################################################################

use strict;
use Getopt::Std;
use File::Find;
use FileHandle::Unget;

# command line syntax
my $ProgFile = ( split('/',__FILE__) )[-1];
sub usage {
	print "\n";
	print <<END_USAGE;
PURPOSE
 in Trinity data directory (genome-guided assembly), find read list files
 associated with final contigs

COMMAND LINE SYNTAX
 $ProgFile [-h]
 $ProgFile trinity-dir > contig2readlist.tsv

arguments:
 trinity-dir  sequence file in fastq format
 contig2readlist.tsv
              translation table with contig ID -> readlist file

options:
 -h           show this command-line syntax description and exit
END_USAGE
	print "\n";
	exit 1;
}

# command-line interface
our %mopt;
getopts('hx:',\%mopt) or &usage();
if ($mopt{h} or int(@ARGV)<1) { &usage() }
our ($tpath, $fseq)=@ARGV;


# recourse on directory content : primary contig set
our %dicrlst;
sub find_do {
	if (m/Trinity\.fasta$/) {
		my $freads=$File::Find::name;
		$freads=~s/\.out\.Trinity\.fasta//; # file suffix with .out.Trinity.fasta
		$freads=~s/\.out\/Trinity\.fasta//; # Trinity.fasta in .out/
		foreach (&read_headers($_)) {
			if (exists($dicrlst{$_})) {
				warn "WARNING: unspecific contig header information, $_\n  1st $dicrlst{$_}\n  2nd $freads\n";
			} else {
				$dicrlst{$_} = $freads;
			}
		}
	}
}
find(\&find_do,$tpath);

# parse final contig set
if (not defined $fseq) {
	$fseq="$tpath/Trinity-GG.fasta";
}
open(INSEQ,$fseq) or die "ERROR: cannot read file $fseq\n";
my $hSeq=FileHandle::Unget->new(\*INSEQ,'r');
while ($_=&read_seq($hSeq)) {
	my $ctg=$_->{id};
	my $c0phrase=$_->{id}.' '.$_->{descr};
	$c0phrase=~s/^.*?_(c\d)/$1/;

	# merge with primary contig set
	if (exists($dicrlst{$c0phrase})) {
		printf "%s\t%s\n", $ctg, $dicrlst{$c0phrase};
	} else {
		warn "WARNING: missing contig in primary contig set: $_->{id}, $c0phrase\n";
	}
}
close(INSEQ);


################################################################################
# functions

# read sequence entries from fasta file and return list of headers
# arg1: sequence file
# ret:  array of header strings (id+descr)
sub read_headers {
  my ($fseq)=@_;#
  my @header;
  open(IN,$fseq) or die "ERROR: cannot read file $fseq\n";
  while (<IN>) {
    if (m/^>(.+)/) { push @header,$1 }
  }
  close(IN);
  return @header;
}

# read single sequence entry from fasta file handle
# arg1: file handle for sequence input, object FileHandle::Unget
#       FileHandle object mediates buffered reading from non-physical handle
#       (pipe)
# ret:  hash reference to {id=>$id,descr=>$descr,seq=>$sSeq}
sub read_seq {
  my ($hInSeq)=@_;
  while (<$hInSeq>) {
    # initialize sequence entry with fasta header line
    if (m/^>(\S+)(?:\s+(.+))?/) {
      my $pSeq={id=>$1,descr=>$2||''};
      my @seq;
      my $lineo;
      while (defined($_=$lineo=<$hInSeq>) and !m/^>/) {
        chomp;
        push @seq,$_;
      }
      # restore filehandle reading position
      if (length($lineo||'')){ $hInSeq->ungets($lineo) }
      # finish sequence entry
      $pSeq->{seq} = join('',@seq);
      return $pSeq;
    }
  }
  return undef;  # no entry anymore
}
