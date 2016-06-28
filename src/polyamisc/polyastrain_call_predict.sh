#!/usr/bin/env bash
mydir=$(dirname $0)  # program directory
if [ "$mydir" == "${mydir#/}" ] ; then mydir=$(pwd)/$mydir ; fi
inf=$1
tmp=/tmp/polyaspwm$$

function fmtpwmout {
  flen=$1 ; shift
  seqct -t $flen  \
  | perl -e '
    our %dics;
    open($hin,shift());
    while (<$hin>){
      if (m/^#/) { next }
      my @v=split;
      $dics{$v[0]}=$v[1];
    }
    close($hin);
    open($hin,shift());
    our ($s, $m,$mp,$msc);
    while (<$hin>){
      my @v=split;
      chomp $v[-1];
      if ($s ne $v[0]) {
        print  "###\n";
        printf "# sequence %s\n", $s=$v[0];
        printf "# seq_length %d\n", $dics{$s};
        $s=$v[0];
      }
      elsif (length($v[2])) {
        if ($v[5] eq $m and $v[2]==$mp+1) { $mp=$v[2]; next; }
        if ($v[2]==$mp+1 and $v[3]<$msc) { $mp=$v[2]; $msc=$v[3]; $m=$v[5]; next; }
      }
      print;
      $m=$v[5]; $mp=$v[2]; $msc=$v[3];
    }' - $@
}

# work, format
echo '### polya_pwm'
echo '###' pwd=$(pwd)
echo '###' mydir=$mydir
echo '###' inf=$inf
echo '###' calltrim=$mydir/fasta_trimpolya.pl
echo '###' callpwm=$(which PWMscanNt)
echo '###' callseqct=$(which seqct)
echo '#'
cat $inf | tee $tmp.in  \
  | $mydir/fasta_trimpolya.pl $inf  \
  | PWMscanNt $mydir/polyastrain_pwm3.pwm  \
  > $tmp
# format, merge sequence lengths
fmtpwmout $tmp.in $tmp

# tidy up
rm ${tmp}*
