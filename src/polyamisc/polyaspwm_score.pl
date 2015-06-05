#!/usr/local/bin/perl
use strict; #use warnings;  # OK 19700101
use Getopt::Std;

our $wsize=15;
our $pk_shift=22;  # ad hoc 25
our $sc_thresh=0.05;
our $colp=2;
our $colsc=3;

# command line usage description
my $ProgFile = ( split('/',__FILE__) )[-1];
sub usage {
  print  "\n";
  print  <<END_USAGE;
PURPOSE
 produce score plot from output of polyaspwm_predict.pl

COMMAND LINE SYNTAX
 $ProgFile [-h]
 $ProgFile [-v] polyaspwm.csv

arguments:
 *.csv      predictiction output from polyaspwm_predict.pl
            (or polyastrain_call_predict.sh)

options:
 -c n       table column that holds the score values, default $colsc (computational)
 -h         print command line usage description and exit
END_USAGE
  print  "\n";
  exit 1;
}

# command-line interface
our %mopt;
getopts('c:h',\%mopt) or &usage();
if ($mopt{h} or int(@ARGV)<1) { &usage() }
if ($mopt{c}) { $colsc=$mopt{c} }
our ($fseq)=@ARGV;


open(IN,$fseq) or die "ERROR: cannot read file $fseq";
our $s=undef;
our $p-1;
our @cand;
while (<IN>) {
  if (m/^#/){ undef $s; print; next; }
  chomp;
  my @v=split/\t/;
  # close previous sequence
  if (length($s) and $v[0] ne $s){
    foreach ( ; int(@cand); ++$p){ &score_eval($p) }
    $p=-1;
    @cand=();
  }
  # next item in current sequence
  $s=$v[0];
  if (int(@cand)) {
    my $pn=&min($v[$colp]+$pk_shift-2*$wsize,$cand[-1][$colp]+$pk_shift+2*$wsize);
    push @cand, \@v;
    foreach ( ; $p<$pn; ++$p){ &score_eval($p) }
  }
  else{
    $p=&max(($p>0)?$p:(),$v[$colp]+$pk_shift-2*$wsize);
    my $pn=$v[$colp]+$pk_shift-2*$wsize;
    push @cand, \@v;
    foreach ( ; $p<$pn; ++$p){ &score_eval($p) }
  }
}

foreach ( ; int(@cand); ++$p){ &score_eval($p) }

exit;


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


sub gaussval {
  my ($my,$sig,$x) = @_;
  if ($sig<=0) { return undef }
  # calculate value
  my $g = 0.3989422804 * exp -(($x-$my)**2 / (2 * $sig**2));
  return $g;
}


sub score_eval {
  my ($p) = @_;
  while (int(@main::cand) and $main::cand[0][$main::colp]+$pk_shift+2*$main::wsize<$p) { shift @main::cand }
  my $sc=0;
  foreach my $pCand (@main::cand) {
    $sc += &gaussval($pCand->[$main::colp]+$main::pk_shift,$main::wsize,$p)
       * 2**$pCand->[$main::colsc];
  }
  if ($sc>$main::sc_thresh){
    printf "%s\t%d\t%.3f\t%s\n", $s, $p, $sc,
      join(',',map{ $_->[$main::colp] } @main::cand);
  }
}
