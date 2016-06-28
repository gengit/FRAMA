#!/usr/bin/env bash
mydir=$(dirname $0)  # program directory
inf=$1
tmp=/tmp/polyaspwmimg$$

function pwmout_sect1st {
  perl -ne '
    BEGIN{ our $sct }
    if (m/^# sequence /) { if($sct){ last } ++$sct }
    print;
    ' $@
}

function pwmout_len1st {
  perl -ne 'if (m/^# seq_length (\d+)/) { print $1,"\n"; exit; }' $@
}

# work, format
l=$(cat $inf | tee $tmp | pwmout_len1st)
echo mydir=$mydir
echo tmp=$tmp
echo l=$l
pwmout_sect1st $tmp  \
| $mydir/polyaspwm_score.pl -  \
| awk '{ printf "%s\t%s\n",$2,$3 }'  \
| Plot.pl -plot -ReprType=column -rangex=0..$l -outstump=polyaspwm -

# tidy up
rm $tmp
