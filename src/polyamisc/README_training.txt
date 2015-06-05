# package polyamisc
#
# copyright (c)
#   FLI Jena, Genome Analysis Group, 2014
# author
#   K Szafranski, szafrans@fli-leibniz.de
#
# notes for training work
#

# insert true/false - only for training phase
# cmp. analogous code in polyastrain_call_predict.sh
function fmtpwmout {
  perl -ne '
    if(m/^#/){ next }
    if(m/total sequence length (\d+)/){ $l=$1; next; }
    @v=split/\t/;
    if(length($v[2])){
      $v[2]=$l-$v[2];
      splice @v,3,0, ($v[2]<=60)?"1":"0";
      if ($v[6] eq $m and $v[2]==($mp-1)) { $mp=$v[2]; next; }
      $m=$v[6]; $mp=$v[2];
      print join(qq(\t),@v);
    }
    '
}

# prediction work
wdir=tmp/cp-201301101652_mull_goldenref
cdir=code/pwmsuite/bin
stump=tmp/polyaspwm
(
  for f in $wdir/*/ref_*.fa ; do
    seqct $f
    ${stump}_trimpolya.pl $f  \
    | PWMscanNt ${stump}_pwm3.pwm
  done
) | fmtpwmout | tee ${stump}_out.csv | wc -l

# evaluate list, for PWM adjustment
perl -ne '
  chomp;
  @v=split/\t/;
  $dics{$v[0]}=1;
  if ($v[0] ne $s){ $s=$v[0]; foreach (keys %dicm){ $dicmglob{$_}+=$dicm{$_} } %dicm=(); }
  if($v[3]==1){ ++$dic{$v[6]}[0]; $dicm{$v[6]}=1; }
  else{ ++$dic{$v[6]}[1] }
  END{
    foreach (sort keys %dic){
      my $oddr=$dic{$_}[1]?$dic{$_}[0]/$dic{$_}[1]:0;
      my $usg=$dicmglob{$_}/int(keys %dics);
      printf "%s\t%d\t%d\t%.3f\t%.3f\t%.3f\n", $_,@{$dic{$_}},
        $oddr, $usg, log(sqrt($oddr*$usg))+0.404;
    }
  }
  ' ${stump}_out.csv

# evaluate list, global scoring
${stump}_revscore.pl -c4 ${stump}_out3.csv | grep ^ACO2_Graumull  \
| awk '{ printf "%s\t%s\n",$2,$4 }' | Plot.pl -plot -ReprType=column -
display /home/lamarck/szafrans/stdin_plot.png &
