#!/usr/bin/env bash

if [ -s $1\_cl_clusters ]
then
    mydir=$(dirname $1)
    cdbyank $1\.cidx < $1\.singletons > $1\.singletons.fa
    find $mydir -name 'contigs' -exec cat {} + > $mydir/all-contigs.fa
    find $mydir -name 'align' -exec cat {} + > $mydir/all-align.fa
    find $mydir -name 'ACE' -exec cat {} + > $mydir/all-ACE.fa
    cat $1\.singletons.fa $mydir/all-contigs.fa > $2
    perl tgicl.pl $mydir/all-ACE.fa > $2.clstr
    rm -rf $mydir/asm_*
else 
    echo "TGICL: Nothing to do."
    cp $1 $2
fi
