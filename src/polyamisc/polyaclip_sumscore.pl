#!/usr/bin/env perl
################################################################################
#
#  polyas_sumscore.pl
#  score poly(A) evidence from different signal sources
#
#  copyright (c)
#    FLI Jena, Genome Analysis Group, 2014
#  author
#    Karol Szafranski, szafrans@fli-leibniz.de
#
################################################################################

use strict; use warnings;  # OK 20140630
use Getopt::Std;
use FileHandle::Unget;

our $co_pas=-0.05;
### ahpp
our $w_pas=1.18;  # 1.18 0.335 0.285 0.99 
our $w_aaa=0.335;
our $w_cdrop=0.285;
our $w_3end=0.99;

# command line syntax
my $ProgFile = ( split('/',__FILE__) )[-1];
sub usage {
  print "\n";
  print <<END_USAGE;
PURPOSE
 score poly(A) evidence from different signal sources

COMMAND LINE SYNTAX
 $ProgFile [-h]
 $ProgFile [options] pwmpolyas.csv polya.csv covdrop.csv 3homol.csv > out.csv

arguments:
 pwmpolyas.csv
            scoring table of poly(A) signals
 polya.csv  scoring table of poly(A) evidence from reads
 covdrop.csv
            scoring table of coverage drop
 3homol.csv scoring table of 3p homology

options:
 -h         show this command-line syntax description and exit
 -a f       set weight to poly(A) read scoring, default $w_aaa
 -A f       multiply default weight to poly(A) read scoring
 -d f       set weight to coverage drop scoring, default $w_cdrop
 -D f       multiply default weight to coverage drop scoring
 -e f       set weight to homologous end scoring, default $w_3end
 -E f       multiply default weight to homologous end scoring
 -p f       set weight to PAS scoring, default $w_pas
 -P f       multiply default weight to PAS scoring
 -q f       set lower cutoff to PAS scoring, default $co_pas
END_USAGE
  print "\n";
  exit 1;
}

# command-line interface
our %mopt;
getopts('ha:A:d:D:e:E:p:P:q:',\%mopt) or &usage();
if ($mopt{h} or int(@ARGV)<4) { &usage() }
if (exists($mopt{a}) and defined($mopt{a})) { $w_aaa=$mopt{a} }
if (exists($mopt{d}) and defined($mopt{d})) { $w_cdrop=$mopt{d} }
if (exists($mopt{e}) and defined($mopt{e})) { $w_3end=$mopt{e} }
if (exists($mopt{p}) and defined($mopt{p})) { $w_pas=$mopt{p} }
if (exists($mopt{q}) and defined($mopt{q})) { $co_pas=$mopt{q} }
if ($mopt{A}) { $w_aaa*=$mopt{A} }
if ($mopt{D}) { $w_cdrop*=$mopt{D} }
if ($mopt{E}) { $w_3end*=$mopt{E} }
if ($mopt{P}) { $w_pas*=$mopt{P} }


# initialize input interfaces to feature lists
# build dictionary @Fls,%dicFt2Fls -> list of FeatureList's
our (@Fls,%dicFt2Fls);
my ($fpas,$faaa,$fcdrop,$f3end) = @ARGV;
our $pPas = ListPwmpolyas->new($fpas);
if (not $pPas) { warn "WARNING: unable to read pwmpolyas list, $fpas - ignored\n" }
our $pAaa = ListPolya->new($faaa);
if (not $pAaa) { warn "WARNING: unable to read polya evidence list, $faaa - ignored\n" }
our $pCdrop = ListCovdrop->new($fcdrop);
if (not $pCdrop) { warn "WARNING: unable to read coverage-drop evidence list, $fcdrop - ignored\n" }
our $p3end = List3homol->new($f3end);
if (not $p3end) { warn "WARNING: unable to read 3p homology evidence list, $f3end - ignored\n" }
my ($pStrHave) = grep{ $_ } $pPas, $pAaa, $pCdrop, $p3end;
if (not $pStrHave) { die "FATAL: unable to read feature list\n" }
foreach my $ppLs (\$pPas, \$pAaa, \$pCdrop, \$p3end) {
  if (not $$ppLs) {
    $$ppLs = ListFake->new($pStrHave->preview_id());
    if (not $$ppLs) { die "ERROR: unable to install fake feature list\n" }
  }
}


# start output with a header
print  "###\n";
printf "# sequence %s\n", $pPas->preview_id();
print  "#\n";

# loop over list stream
my $p=0;
while (defined($p=&next_pos(++$p))) {

  # weigh feature scores
  # ad-hoc parameters were 1.0 1.2 0.36 1.0
  # ad-hoc(+) parameters are 1.00 0.27 0.33 1.05
  # ad-hoc(++)=opt1 parameters are 1.11 0.32 0.33 1.03
  # martin1 parameters 0.79 1.07 0.53 1.59
  # martin1(+))) parameters 0.79 1.07 0.53 1.59
  # martin2 parameters 1.33 0.22 0.44 0.66
  # opt2 parameters are 1.11 0.35 0.30 1.03 (feature peaks re-calibrated; polyareads normalized)
  my $sc_pas = &max(0,$pPas->score_pos($p)-$co_pas) *$w_pas;
  my $sc_aaa = $pAaa->score_pos($p) *$w_aaa;
  my $sc_cdrop = &max( $pCdrop->score_pos($p) ,0) *$w_cdrop;
  my $sc_3end = $p3end->score_pos($p) *$w_3end;

  # create a combined score
  printf "%s\t%d\t%.3f\t%s\t%s\t%s\t%s\n", $pPas->preview_id()||'undef', $p,
    $sc_pas+$sc_aaa+$sc_cdrop+$sc_3end, $sc_pas, $sc_aaa, $sc_cdrop, $sc_3end;
}


# scrawl position by position in sequence by sequence
# individual lists may lack porisition ranges of a sequence
sub next_pos {
  my ($p) = @_;
  my $pret=undef;
  foreach my $pLs ($pPas, $pAaa, $pCdrop, $p3end) {
    $pLs->preview_pos($p);
    if (defined($pLs->{pos})) {
      if ($pLs->{pos}<$p) {
        printf "WARNING: pos out of order, requested %d, next %s\n", $p, $pLs->{pos}||'undef';
      } else {
        $pret=$p;
      }
    }
  }
  return $pret;
}


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


################################################################################
# interface to input list
#
# data:
#  seqid   current sequence id
#  pos     current sequence position
#  buffer  reference to list of buffered lines
#
package ListMain;
use Exporter qw(import);
our @EXPORT = qw (
  new preview_id preview_pos next_seq
  );

use FileHandle;

sub new {
  my ($pkg,$f) = @_;
  my $this = { path=>$f, hIn=>FileHandle->new($f,'r') };
  $this->{hIn} or die "ERROR: cannot open feature list file $f";
  $this->{seqid} = do{ {
    my $l=$this->{hIn}->getline();
    if ($this->{hIn}->eof()) {
      $this->{seqid}=undef;
      $this->{pos}=undef;
      return undef;
    }
    chomp $l;
    if ($l=~m/^#\s*sequence (\w+)/) { $1 }
    else{ redo }
  } };
  bless $this, $pkg;
}

sub preview_id {
  my ($this) = @_;
  return $this->{seqid};
}

sub preview_pos {
  my ($this,$p) = @_;
  if (!defined($this->{seqid})) { return undef }
  if (exists($this->{seqend}) and $this->{seqend}) { return undef }
  while (
    defined($this->{seqid}) and
    !(exists($this->{seqend}) and $this->{seqend}) and
    (!defined($this->{pos}) or $this->{pos}<$p) and
    my $l=$this->{hIn}->getline()
  ) {
    chomp $l; $this->{buffer}=$l;
    if ($this->{hIn}->eof() or $l=~m/^###/) {
      $this->{seqend}=1;
      $this->{pos}=undef;
      return undef;
    }
    if ($l=~m/^#/) { next }
    $this->{pos} = (split/\t/,$l)[$this->column_pos()];
  }
  return $this->{pos};
}

sub next_seq {
  my ($this) = @_;
  $this->{pos}=undef;
  delete $this->{seqend};
  return $this->{seqid} = do{ {
    my $l=$this->{hIn}->getline();
    if ($this->{hIn}->eof()) {
      $this->{seqid}=undef;
      $this->{pos}=undef;
      return undef;
    }
    chomp $l; $this->{buffer}=$l;
    if ($l=~m/^#\s*sequence (\w+)/) { $1 }
    else{ redo }
  } };
}

sub score_pos {
  my ($this,$p) = @_;
  if (!defined($this->{seqid})) { return undef }
  if (!defined($this->{pos}) or $p<$this->{pos}) {
    return $this->default_score();
  }
  else{
    return (split/\t/,$this->{buffer})[$this->column_score()];
  }
}


################################################################################
# interface to input list
package ListPwmpolyas;
#use Exporter qw(import);
#ListMain->import();

sub new { return ListMain::new(@_); }

sub preview_id { return ListMain::preview_id(@_); }

sub preview_pos { return ListMain::preview_pos(@_); }

sub next_seq { return ListMain::next_seq(@_); }

sub column_pos { 1; }

sub column_score { 2; }

sub default_score { 0.05; }

sub score_pos { return ListMain::score_pos(@_); }


################################################################################
# interface to input list
package ListPolya;
#use Exporter qw(import);
#require ListMain; ListMain->import();

sub new { return ListMain::new(@_); }

sub preview_id { return ListMain::preview_id(@_); }

sub preview_pos { return ListMain::preview_pos(@_); }

sub next_seq { return ListMain::next_seq(@_); }

sub column_pos { 1; }

sub column_score { 2; }

sub default_score { 0.0; }

sub score_pos { return ListMain::score_pos(@_); }


################################################################################
# interface to input list
package ListCovdrop;
#use Exporter qw(import);
#ListMain->import();

sub new { return ListMain::new(@_); }

sub preview_id { return ListMain::preview_id(@_); }

sub preview_pos { return ListMain::preview_pos(@_); }

sub next_seq { return ListMain::next_seq(@_); }

sub column_pos { 1; }

sub column_score { 2; }

sub default_score { 0.0; }

sub score_pos { return ListMain::score_pos(@_); }


################################################################################
# interface to input list
package List3homol;
#use Exporter qw(import);
#ListMain->import();

sub new { return ListMain::new(@_); }

sub preview_id { return ListMain::preview_id(@_); }

sub preview_pos { return ListMain::preview_pos(@_); }

sub next_seq { return ListMain::next_seq(@_); }

sub column_pos { 1; }

sub column_score { 2; }

sub default_score { 0.0; }

sub score_pos { return ListMain::score_pos(@_); }


################################################################################
# interface to input list
package ListFake;
#use Exporter qw(import);
#ListMain->import();

sub new {
  my ($pkg,$s) = @_;
  my $this = { seqid=>$s, seqend=>1 };
  bless $this;
}

sub preview_id { return ListMain::preview_id(@_); }

sub preview_pos { return ListMain::preview_pos(@_); }

sub next_seq { return ListMain::next_seq(@_); }

sub default_score { 0.0; }

sub score_pos { return ListMain::score_pos(@_); }
