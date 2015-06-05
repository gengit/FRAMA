#!/usr/bin/env perl

use strict;
use warnings;

=pod

=head1 INFO

Martin Bens, bensmartin@gmail.com
2014-08-27

=head1 DESCRIPTION

Annotated GENSCAN predictions. 

Hits file must contain:
    QUERY_ID, STRAND, START_IN_TARGET, STOP_IN_TARGET, IDENTITY, COVERED_QUERY, QUERY_SYMBOL

Additionally, use known CDS regions to annotated result in the follow format
    START, END, STRAND, SYMBOL, ACCESSION

=head1 OPTIONS

=over 8

=item B<-genscan>

Genscan result as produced by GENSCAN.

=item B<-hits> 

QueryID, Strand, StartInTarget, EndInTarget, Identity, CoveredQuery, Symbol.

=item B<-cds>

Already known coding regions (start, end, strand, symbol, accession).

=item B<-rev> 

Specify length of contig if blast was performed on reverse complemented
sequence. 

=item B<-head/nohead> 

Header line.

=back

=cut

use Getopt::Long;
use Pod::Usage;
use Bio::Tools::Genscan;
use Data::Dumper;
use Bio::Range;

pod2usage(-message => "\n\tNo arguments\n", -verbose => 1) if (@ARGV == 0);

# cut -f1,8,11,12,16,17 blast.csv > blast_preprocessed.csv
# repCol.pl -table1 blast_preprocessed.csv -table2 $RNA_HS/acc2sym.csv -tkey 1 -tvalue 2 -column 1 -add > blast_preprocessed_symbol.csv
# perl annotated_genscan.pl -genscan CDS_genscan.txt -blast blast_preprocessed_symbol.csv -rev 4820

# Blast
use constant TARGET_ID    => 0;
use constant TARGET_START => 1;
use constant TARGET_END   => 2;
use constant IDENTITY     => 3;
use constant SYMBOL       => 4;

use constant MAX_INT => 9**9**9;

my $man  = 0;
my $help = 0;
my ($genscan_file, $blast_file, $cds_file);
my $debug = 0;
my $rev   = 0;
my $head  = 1;
GetOptions(
    'genscan|g=s' => \$genscan_file,
    'hits|b=s'    => \$blast_file,
    'cds|c=s'     => \$cds_file,
    'debug|d'     => \$debug,
    'rev|r=s'     => \$rev,
    'head!'       => \$head,
    'help|?'      => \$help,
    'man|m'       => \$man
) or pod2usage(1);

pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;

if (!-e $genscan_file || !-e $blast_file) {
    die
      "\n\tSpecify genscan result and blast result (blast result as produced by wublast_average.pl\n";
}

my %cds;
if ($cds_file) {
    unless (-e $cds_file) {
        print STDERR "Could not find CDS file\n";
    }
    print STDERR "Reading CDS regions...\n" if ($debug);

    open my $fh, "<", $cds_file or die $!;
    while (<$fh>) {
        next if (/^#/);
        chomp;
        my @e = split;
        $cds{$e[4]} = {
            range      => Bio::Range->new(-start => $e[0], -end => $e[1], -strand => $e[2]),
            symbol     => $e[3],
            transcript => $e[4],
            annotation => 1
        };
    }
    close $fh;
}

my ($predicted_cds, $prom_region) = getGenscanResults($genscan_file);

# overlap with existing CDS?
my %found_annotated_cds;
my @final_cds_regions;
my $assignSymbol = 0;
for my $prediction (@$predicted_cds) {
    my $overlap = 0;
    for my $annotation (values %cds) {
        my ($start, $end, $strand) = $prediction->intersection($annotation->{range});
        if ($start) {

            # might have been shown an overlap previously (multiple prediction for on CDS)
            if (!exists $found_annotated_cds{$annotation->{transcript}}) {
                $found_annotated_cds{$annotation->{transcript}} = 1;
                push @final_cds_regions, $annotation;
            }
            $overlap = 1;
        }
    }
    unless ($overlap) {
        $assignSymbol = 1;    # at least one unknown CDS to annotat
        push @final_cds_regions, {range => $prediction};
    }
}

# add remaining CDS without overlap with predicted genes
for my $annotation (values %cds) {
    push @final_cds_regions, $annotation
      unless (exists $found_annotated_cds{$annotation->{transcript}});
}

my (@genes, @transcripts);
if ($assignSymbol) {
    my @possible_cds = sort { $a->{range}->start <=> $b->{range}->start } @final_cds_regions;

    # blast: accession, start, end, identity, coverage, symbol
    my @blast = @{getBlastHash($blast_file)};
    if ($rev) {
        @blast = map {
            [
                $_->[TARGET_ID], $rev - $_->[TARGET_END], $rev - $_->[TARGET_START],
                $_->[IDENTITY],  $_->[SYMBOL]
            ]
        } @blast;
    }
    @blast = sort { $a->[TARGET_START] <=> $b->[TARGET_START] } @blast;

    # find best blast hit which overlaps predicted CDS
    # if no such blast hit exists, we assume wrong prediction
    my %best_cds_per_gene;
    my $counter = 0;
  NEXT_CDS: for my $c_cds (@possible_cds) {

        if (exists $c_cds->{symbol}) {

            # already known
            $best_cds_per_gene{$c_cds->{symbol}} = {
                length => MAX_INT,
                cds    => $counter,
                acc    => $c_cds->{transcript},

            };
        } else {

            # look for blast overlap

            for my $c_blast (@blast) {

                # no further overlapping blast results (hit start > cds end)
                if ($c_blast->[TARGET_START] > $c_cds->{range}->end) {
                    $counter++;
                    next NEXT_CDS;
                }

                # overlap between blast hit and genscan prediction
                my ($start, $stop, undef) = $c_cds->{range}->intersection(
                    Bio::Range->new(
                        -start => $c_blast->[TARGET_START],
                        -end   => $c_blast->[TARGET_END]
                    )
                );

                # for each gene, keep cds with longest overlap (or highest identity if equal length)
                if ($start) {
                    my $length = $stop - $start;

                    my $gene_symbol = $c_blast->[SYMBOL];
                    if (exists $best_cds_per_gene{$gene_symbol}) {
                        if ($best_cds_per_gene{$gene_symbol}->{length} < $length) {
                            $best_cds_per_gene{$gene_symbol} = {
                                length   => $length,
                                identity => $c_blast->[IDENTITY],
                                cds      => $counter,
                                acc      => $c_blast->[TARGET_ID]
                            };
                        } elsif ($best_cds_per_gene{$gene_symbol}->{length} == $length
                            && $best_cds_per_gene{$gene_symbol}->{identity} > $c_blast->[IDENTITY])
                        {
                            $best_cds_per_gene{$gene_symbol} = {
                                length   => $length,
                                identity => $c_blast->[IDENTITY],
                                cds      => $counter,
                                acc      => $c_blast->[TARGET_ID]
                            };
                        }
                    } else {
                        $best_cds_per_gene{$gene_symbol} = {
                            length   => $length,
                            identity => $c_blast->[IDENTITY],
                            cds      => $counter,
                            acc      => $c_blast->[TARGET_ID]
                        };
                    }
                }
            }
        }
        $counter++;
    }

    # sort by length and
    for my $key (
        sort { $best_cds_per_gene{$b}->{length} <=> $best_cds_per_gene{$a}->{length} }
        keys %best_cds_per_gene
      )
    {
        my $value = $best_cds_per_gene{$key};
        next if (exists $final_cds_regions[$value->{cds}]->{symbol});
        $final_cds_regions[$value->{cds}]->{symbol}     = $key;
        $final_cds_regions[$value->{cds}]->{transcript} = $value->{acc};
        $final_cds_regions[$value->{cds}]->{annotation} = $value->{annotation};
    }
}

if ($head) {
    print join "\t", "#feature", "symbol", "accession", "start", "end", "strand", "annotation";
    print "\n";
}
for (my $i = 0; $i < @final_cds_regions; $i++) {
    my $predicted = $final_cds_regions[$i]->{annotation} ? 1 : 0;

    unless ($final_cds_regions[$i]->{symbol}) {
        $final_cds_regions[$i]->{symbol}     = "UNKNOWN";
        $final_cds_regions[$i]->{transcript} = "UNKNOWN";
    }

    print join "\t", "CDS", $final_cds_regions[$i]->{symbol}, $final_cds_regions[$i]->{transcript},
      $final_cds_regions[$i]->{range}->start,  $final_cds_regions[$i]->{range}->end,
      $final_cds_regions[$i]->{range}->strand, $predicted;
    print "\n";
    if ($prom_region->[$i]) {
        print join "\t", "PRO", $final_cds_regions[$i]->{symbol},
          $final_cds_regions[$i]->{transcript}, $prom_region->[$i]->start, $prom_region->[$i]->end,
          $final_cds_regions[$i]->{range}->strand, 0;
        print "\n";
    }
}

sub getGenscanResults {
    my ($out_genscan_file) = @_;
    my $genscan = Bio::Tools::Genscan->new(-file => $out_genscan_file);

    my @predictions;
    eval {
        while (my $gene = $genscan->next_prediction()) {
            push @predictions, $gene;
        }
    };
    if ($@) {
        print STDERR "Could not retrieve genscan result. Skipping genscan!\n";
        return;
    }

    my @cds_regions;
    my @prom_region;
    for my $gene (@predictions) {
        next unless ($gene->exons);
        my (@sorted) = sort { $a->location->start <=> $b->location->start } $gene->exons;
        my $start = $sorted[0];
        my $end = $sorted[$#sorted];

        push @cds_regions,
          Bio::Range->new(
            -start  => $start->location->start,
            -end    => $end->location->end,
            -strand => $gene->strand
          );

        my $range;
        if ($gene->promoters) {
            my @sorted_exon = sort { $a->start <=> $b->start } $gene->promoters();
            if ($gene->strand == -1) {

                # most 5' promotor (in contig orienation)
                $start = $sorted_exon[0]->start;
                $end   = $sorted_exon[0]->end;
                $range = Bio::Range->new(-start => $start, -end => $end, -strand => $gene->strand);
            } else {

                # most 3' promotor (in contig orienation)
                $start = $sorted_exon[$#sorted_exon]->start;
                $end   = $sorted_exon[$#sorted_exon]->end;
                $range = Bio::Range->new(-start => $start, -end => $end, -strand => $gene->strand);
            }
        }
        push @prom_region, $range;

    }
    $genscan->close();

    return \@cds_regions, \@prom_region;
}

sub getBlastHash {
    my ($file) = @_;
    my @result;
    open my $fh, "<", $file or die $!;
    while (<$fh>) {
        next if (/^#/);
        my @e = split;
        push @result, [$e[TARGET_ID], $e[TARGET_START], $e[TARGET_END], $e[IDENTITY], $e[SYMBOL]];
    }
    close $fh;

    return \@result;
}

1;

__END__

