
ssh gen102

cdir=/gen/fly/biosw/polyamisc
fctg=/misc/vulpix/data/mbens/projects/trinity/2014-03-14_Heterocephalus_glaber_v3/output/sequences-mRNA.fa
fdout=/home/lamarck/mbens/vulpix/projects/2014-06-23_polyA/features/trinity_artificial/

# gene in datasets
gene_choose CIZ1
# effectively:
#g=CIZ1
ctg=$(gettrinity_spec_accession $g | sort +1 | head -1 | cut -f1)
$ctg | fasta_slcid.pl - $fctg > tmp/polya_wk/${g}_seq.fa

# score poly(A) signal
$cdir/polyaspwm_predict.pl tmp/polya_wk/${g}_seq.fa | $cdir/polyaspwm_score.pl -  \
> tmp/polya_wk/${g}_polyas.csv

# score poly(A) events found in readings
cdir=/gen/fly/biosw/polyamisc
~mbens/dev/pipeline/src/bam_polyA.pl -i $fdout/$ctg.bam  \
| $cdir/polyaread_score.pl -  \
> tmp/polya_wk/${g}_polya.csv

# raw read coverage
cdir=/gen/fly/biosw/polyamisc
if [ -s "tmp/polya_wk/${g}_ocover.csv" ]
  then callcover="cat tmp/polya_wk/${g}_ocover.csv"
  else callcover="bamtools coverage -in $fdout/$ctg.bam | tee tmp/polya_wk/${g}_ocover.csv"
fi
eval $callcover  \
| $cdir/polyacov_coverage.pl -  \
> tmp/polya_wk/${g}_cover.csv
eval $callcover  \
| $cdir/polyacov_dropcov_simple.pl -  \
> tmp/polya_wk/${g}_dropcov.csv

# score homology matches to 3p end
$cdir/csv_slccontig.pl $ctg $fdout/summary_human_aln.csv  \
| $cdir/polyahomol_score.pl -  \
> tmp/polya_wk/${g}_3homol.csv

# global score
$cdir/polyaclip_sumscore.pl tmp/polya_wk/${g}_polyas.csv tmp/polya_wk/${g}_polya.csv tmp/polya_wk/${g}_dropcov.csv tmp/polya_wk/${g}_3homol.csv | more  \
> tmp/polya_wk/${g}_polyasum.csv
