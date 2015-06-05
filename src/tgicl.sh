#!/usr/bin/env bash

if [ -s $2\_cl_clusters ]
then
  cdbyank $1\.cidx < $2\.singletons > $2\.singletons.fa
  find $3 -name 'contigs' -exec cat {} + > $3/all-contigs.fa
  find $3 -name 'align' -exec cat {} + > $3/all-align.fa
  find $3 -name 'ACE' -exec cat {} + > $3/all-ACE.fa
  cat $2\.singletons.fa $3/all-contigs.fa > $4
  mv $3/Trinity.tmp.fa_cl_clusters $4.clstr
  rm -rf $3/asm_*
else 
  echo "TGICL: Nothing to do."
  cp $1 $4
fi
