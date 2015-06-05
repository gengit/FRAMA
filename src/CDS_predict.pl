#!/usr/bin/env perl

=pod

=head1 INFO
Martin Bens, bensmartin@gmail.com
2014-07-22

=head1 DESCRIPTION

=head1 OPTIONS

=over 8

=item B<-contig> contig.fa

Sequence with unknown CDS.

=item B<-ortholog> orth.fa

Ortholog.

=item B<-seleno> 0

Number of selenopositions [default: 0]

=item B<-cds> "1 10"

Start and end of CDS (1-based). End must not include Stop-Codon.

=item B<-species> Hsapiens

Abbreviated species name (e.g. Homo sapiens -> Hsapiens) of provided ortholog.

=item B<-msa> aln.aln

Multiple species alignment. Mark CDS sequences which should be used for
inference with prefix "cds:". Furthermore, an abbreviated species name (Homo
sapiens = Hsapiens) is necessary. A valid sequence name looks like
"cds:Hsapiens:NM_000014"

=item B<-msa-format> fasta

Format of -msa (fasta, clustalw,...)

=item B<-unaligned> 

Specify if -msa contains only unaligned transcripts (in fasta format).

=item B<-intron-length> 50

Minimum required length of gap in orthologous cds to examine for splice site
(GT..AG)

=item B<-predictions>

Ignore RefSeq-Acc starting with "XM" for CDS inference.

=item B<-out-aln> out.aln

Output of alignment including contig and ortholog.

=item B<-out-genscan> genscan.txt

Output of genscan result if computed.

=item B<-all> 

Force to calculate CDS prediction with all methods (otherwise stops after first
successfull prediction)

=item B<-genscan-matrix>



=back


=cut

use strict;
use warnings;

use Getopt::Long;
use File::Basename;
use Data::Dumper;
use Pod::Usage;
use File::Spec::Functions qw(catfile catdir);
use File::Path qw(make_path remove_tree);
use File::Temp qw(tempfile);
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::DB::Fasta;
use Bio::Tools::Genscan;

use FindBin;
use lib "$FindBin::Bin/../lib/perl";
use StringHelper;

use Bio::Tools::CodonTable;

use constant MAFFT_LENGTH => 30000;

my %translateStatus = (
    "001" => "partial;intron",
    "101" => "partial_3prime;intron",
    "011" => "partial_5prime;intron",
    "111" => "full_length;intron",
    "000" => "partial",
    "100" => "partial_3prime",
    "010" => "partial_5prime",
    "110" => "full_length",
    "002" => "partial;pseudogene",
    "102" => "partial_3prime;pseudogene",
    "012" => "partial_5prime;pseudogene",
    "112" => "full_length;pseudogene",
);

my $codon_table = Bio::Tools::CodonTable->new(-id => 1);

my %opt = (
    predictions => 0,
    reindex     => 0,

    'msa-format'      => 'fasta',
    'predictions'     => 0,
    'all'             => 0,
    'intron-length'   => 50,
    'selenopositions' => 0,

    'help'      => 0,
    'man'       => 0,
    'unaligned' => 0,
);

pod2usage(-message => "\n\tNo arguments. See -help.\n") if (@ARGV == 0);

GetOptions(
    \%opt,             "help|h",         "man",           "contig=s",
    "ortholog=s",      "seleno=i",       "msa=s",         "msa-format=s",
    'intron-length=i', 'predictions',    'cds=s',         'out-aln=s',
    'out-cds=s',       'out-prot=s',     'out-genscan=s', 'all',
    'species=s',
    'unaligned',       'help',           'compareTo=s', 'genscan-matrix=s'
) or pod2usage(1);

pod2usage(1) if $opt{help};

# mandatory
unless ($opt{contig}) {
    pod2usage(-message => "\n\tNot enough arguments. See -help.\n");
}

#unless (-e $opt{contig} && -e $opt{ortholog}) {
#    die "\n\tContig or ortholog file not found.\n\n";
#}

#my ($cds_start, $cds_end) = split " ", $opt{cds};

#unless ($cds_start && $cds_end) {
#    die "\n\t-cds should contain somethine like \"20 50\"\n\n";
#}

my $io = Bio::SeqIO->new(-file => $opt{contig}, -format => "fasta");
my $contig_seq = $io->next_seq if ($io);

unless ($contig_seq) {
    die "\n\tCould not find contig sequence: $opt{contig}\n\n";
}

#$io = Bio::SeqIO->new(-file => $opt{ortholog}, -format => "fasta");
#my $orth_seq = $io->next_seq if ($io);

#unless ($orth_seq) {
#    die "\n\tCould not find reference sequence: $opt{ortholog}\n\n";
#}

my $contig_id = $contig_seq->id;

#my $orth_id   = $orth_seq->id;

# Header
print "###\n";
print "# sequences " . $contig_seq->id . "\n";
print "# selenopositions " . $opt{seleno} . "\n";
print "##\n\n";

my $aln          = Bio::SimpleAlign->new();
my $add_ortholog = 1;
my @lengths;

my ($fh_align_sequences, $align_sequences) = tempfile(UNLINK => 1);
my $toAlign = Bio::SeqIO->new(
    -file   => ">" . $align_sequences,
    -format => "fasta"
);

my $compareTo_seq;

# get multiple species alignment
if ($opt{msa}) {
    my $failed = 0;
    if (-e $opt{msa}) {
        if ($opt{unaligned}) {
            my $io = Bio::SeqIO->new(-file => $opt{msa}, -format => $opt{'msa-format'});
            while (my $seq = $io->next_seq) {
                if (!$opt{predictions} && $seq->id =~ /XM/) {
                    next;
                }
                if ($seq->id =~ /cds:(.+?):$opt{compareTo}$/) {
                    $compareTo_seq = $seq->clone;
                }

                #                $add_ortholog = 0 if ($seq->id =~ /$orth_id$/);
                $toAlign->write_seq($seq);
                push @lengths, $seq->length;
            }
        } else {

            my $io = Bio::AlignIO->new(
                -file   => $opt{msa},
                -format => $opt{'msa-format'}
            );
            my $tmp_aln = $io->next_aln() if ($io);
            if ($tmp_aln) {

                # CDS of reference transcript already present?
                for ($tmp_aln->each_seq) {
                    if (!$opt{predictions} && $_->id =~ /XM/) {
                        next;
                    }
                    if ($_->id =~ /cds:(.+?):$opt{compareTo}$/) {
                        $compareTo_seq = $_->clone;
                    }

                    #$add_ortholog = 0 if ($_->id =~ /$orth_id/);
                    $aln->add_seq($_);
                }
                push @lengths, $aln->length;
            } else {
                $failed = 1;
            }
        }
    }
    if ($failed) {
        print STDERR "Could not find multiple species alignment: $opt{msa}, $opt{'msa-format'}\n";
        print STDERR "Using reference sequence only!";
    }
}

# add full length reference transcript in any case
#my $aln_orth_seq = Bio::Seq->new(
#    -id  => "full:$opt{species}:$orth_id",
#    -seq => $orth_seq->seq
#);
#$toAlign->write_seq($aln_orth_seq);
#push @lengths, $aln_orth_seq->length;

#my $orth_cds_seq = Bio::Seq->new(
#    -id  => "cds:$opt{species}:$orth_id",
#    -seq => $orth_seq->subseq($cds_start, $cds_end)
#);
#if ($add_ortholog) {
#    $toAlign->write_seq($orth_cds_seq);
#}

# add contig
$toAlign->write_seq($contig_seq);
push @lengths, $contig_seq->length;

my $output_file;
my $out_fh;
if ($opt{'out-aln'}) {
    $output_file = $opt{'out-aln'};
} else {
    ($out_fh, $output_file) = tempfile(UNLINK => 1);
}

unless ($aln) {

    # align all sequences

    my $command;

    # more accurate settings but not suitable for long alignments
    if (max(@lengths) < MAFFT_LENGTH) {
        $command =
          "mafft --localpair --maxiterate 1000 --quiet $align_sequences > $output_file 2> $output_file.err";
    } else {
        $command = "mafft --quiet $align_sequences > $output_file 2> $output_file.err";
    }
    system($command) == 0 || die "\nFailed to compute alignment with MAFFT.\n\n";
    $io = Bio::AlignIO->new(-file => $output_file, -format => "fasta");
    $aln = $io->next_aln if ($io);
} else {

    # preserve existing alignment and add reference, reference cds, contig

    my ($fh1, $alignment_file) = tempfile(UNLINK => 1);
    my $io = Bio::AlignIO->new(
        -file   => ">" . $alignment_file,
        -format => "fasta"
    );
    $io->write_aln($aln);

    my $command;

    # more accurate settings but not suitable for long alignments
    if (max(@lengths) < MAFFT_LENGTH) {
        $command =
          "mafft --localpair --maxiterate 1000 --quiet --add $align_sequences $alignment_file > $output_file 2> /dev/null";
    } else {
        $command =
          "mafft --quiet --add $align_sequences $alignment_file > $output_file 2> /dev/null";
    }
    system($command) == 0 || die "\nFailed adding sequences to existing alignment\n\n";
    $io = Bio::AlignIO->new(-file => $output_file, -format => "fasta");
    $aln = $io->next_aln;

    close $fh1;
}
close $out_fh unless ($opt{'out-aln'});
close $fh_align_sequences;

my $cds_predicted;
if ($aln) {
    $cds_predicted = getCDSbyMSA($aln, $opt{seleno}, $contig_id);
    if (not defined $cds_predicted) {
        warn("\nCDS prediction from alignment failed.\n");
    }
}

if (!defined $cds_predicted || !$cds_predicted || $opt{'all'}) {
    $cds_predicted = getCDSbyGenscan();
    if (not defined $cds_predicted) {
        warn("\nCDS prediction by GENSCAN failed.\n");
    }
}

if (!defined $cds_predicted || !$cds_predicted || $opt{'all'}) {
    $cds_predicted = getCDSbyLongestOrf();
    if (not defined $cds_predicted) {
        warn("\nCDS prediction by EMBOSS failed.\n");
    }
}

# this should not happen...
if (!defined $cds_predicted || !$cds_predicted) {
    die "\nContig $contig_id has no ORF\n\n";
}

sub getIntronString {
    my ($introns) = @_;

    my @out;
    for (@$introns) {
        push @out, $_->[0] . "-" . $_->[1];
    }
    return join "#", @out if (@out);
    return 0;
}

sub validateSequence {
    my ($prot, $possible_orf, $selenopositions) = @_;

    my @selenos;
    if ($selenopositions > 0) {

        # selenoproteins: allow same number of TGA codons as ortholog
        my $seleno_in_contig = 0;
        while ($prot =~ /\*/g) {
            if (scalar @selenos > $selenopositions) {
                return (0, undef);
            }

            if (uc(substr($possible_orf, ($-[0] * 3), 3)) eq "TGA") {
                push @selenos, $-[0];
            } else {
                return (0, undef);
            }
        }
    } else {
        return (0, undef) if ($prot =~ /\*/);
    }

    if (@selenos > 0) {
        return 1, \@selenos;
    } else {
        return 1, undef;
    }
}

sub getCDSbyMSA {
    my ($aln, $selenopositions, $contig_id) = @_;

    unless ($selenopositions) {
        $selenopositions = 0;
    }

    unless ($aln) {
        return undef;
    }

    my ($aligned_contig_seq) =
      grep { $_->id =~ /$contig_id/ } $aln->each_seq();
    unless ($aligned_contig_seq->seq) {
        warn("Trinity contig not found in alignment: $contig_id");
        return undef;
    }

    (my $seq_string = $aligned_contig_seq->seq) =~ s/-//g;

    my %visited;
    for my $current_seq ($aln->each_seq()) {
        next unless ($current_seq->id =~ /^cds:/);

        my (undef, $species, $orth_id) = split ":", $current_seq->display_id;

        # sequence start and end in alignment
        my ($start, $right) = StringHelper::getOverhang($current_seq->seq, "-");
        my $end = ($aln->length - $right - 1);

        # use first base in contig as start (partial 5'prime)
        my $first_codon = uc(substr($aligned_contig_seq->seq, $start, 3));
        if ($first_codon =~ /^-/) {
            my $cut_start = uc(substr($aligned_contig_seq->seq, $start));
            $start += $-[0] if ($cut_start =~ /[ACGT]/);
        }

        # find gaps > 20 with splicing motifs (GT..AG)
        my $introns =
          findIntrons($current_seq->seq, $aligned_contig_seq->seq, $opt{'intron-length'});
        my $intron_string = getIntronString($introns);

        # check all frames
        for my $frame (0 .. 2) {

            my $current_start = $start + $frame;

            # positions in contig
            my $seq_start = $aligned_contig_seq->location_from_column($current_start + 1);
            my $seq_end   = $aligned_contig_seq->location_from_column($end + 1);
            $seq_start = $seq_start ? $seq_start->start   : 1;
            $seq_end   = $seq_end   ? $seq_end->start : 1;

            my $visit_string = $seq_start . "_" . $seq_end . "_" . $intron_string;

            # CDS of ortholog not aligned
            if ($seq_start == $seq_end) {
                $visited{$visit_string} = {
                    ids    => ["$species:$orth_id"],
                    prot   => "",
                    start  => $seq_start,
                    end    => $seq_end,
                    seleno => "NA",
                    exons  => undef,
                    status => undef,
                };
                next;
            }

            if (exists $visited{$visit_string}) {
                if ($visited{$visit_string}) {
                    push @{$visited{$visit_string}->{ids}}, "$species:$orth_id";
                }
                next;
            }

            # 5prime full, 3prime full, no intron
            my $status = [1, 1, 0];

            my $seq_length = $seq_end - $seq_start + 1;
            $seq_length -= $seq_length % 3;

            my $possible_orf = substr($seq_string, $seq_start - 1, $seq_length);
            my $prot = $codon_table->translate($possible_orf);

            # protein must have "M" as first aminoacid in order to have a
            # complete 5'

            $status->[0] = 0 if ($prot !~ /^M/);

            # in frame stop codons?
            my ($valid_cds, $seleno_in_contig) =
              validateSequence($prot, $possible_orf, $selenopositions);

            my $exons;
            if (!$valid_cds && @$introns > 0) {
                (my $intronless_string, $exons) = getIntronLess($aligned_contig_seq, $introns, $current_start, $end);
                $prot = $codon_table->translate($intronless_string);
                ($valid_cds, $seleno_in_contig) = validateSequence($prot, $possible_orf, $selenopositions);
                $status->[2] = 1;
            }

            unless ($seleno_in_contig) {
                $seleno_in_contig = "NA";
            }

            if (!$valid_cds) {

                # pseudogene, assembly error, annotation error, ...
                $visited{$visit_string} = {
                    ids    => ["$species:$orth_id"],
                    prot   => $prot,
                    start  => $seq_start,
                    end    => $seq_end,
                    seleno => $seleno_in_contig,
                    exons  => undef,
                    status => undef,
                };
                next;
            }

            # valid CDS, look for stop codon
            my $next_start;
            if ($exons) {
                $next_start = $exons->[@{$exons}-1]->[1] + 1;
            } else {
                $next_start = $seq_start + $seq_length ;
            }
            my $prot_end   = "";
            unless ($next_start >= length($seq_string)) {
                $prot_end = $codon_table->translate(substr($seq_string, $next_start - 1 ));
            }

            if ($prot_end =~ /\*/) {
                if ($-[0] == 0) {

                    # same end
                    $seq_end = $next_start + 2;
                } else {

                    # diff end, but stop codon => full 3prime
                    $prot .= substr($prot_end, 0, $-[0]);
                    $seq_end = $next_start + 2 + (3 * $-[0]);
                }
            } else {

                # no stop codon, last complete codon => partial 3prime
                $seq_end = length($seq_string);
                $status->[1] = 0;
                $prot .= $prot_end;
            }

            # adjust coordinates of last exons too
            $exons->[@{$exons} - 1]->[1] = $seq_end if ($exons);

            # incomplete start
            # => maybe result of unclean alignemnt
            # => look for start codon upstream
            if ($seq_start > 2) {
                my $start             = ($seq_start - 1) % 3;
                my $seq_before_start  = substr $seq_string, $start, $seq_start - $start - 1;
                my $prot_before_start = $codon_table->translate($seq_before_start);

                # find leftmost "M" without subsequent '*'
                if ($prot_before_start =~ /(M((?:(?!\*).)*)$)/) {
                    $seq_start   = ($-[0] * 3) + $start + 1;
                    $prot        = $1 . $prot;
                    $status->[0] = 1;
                }
            }

            $visited{$visit_string} = {
                start  => $seq_start,
                end    => $seq_end,
                status => $status,
                ids    => ["$species:$orth_id"],
                exons  => $exons,
                seleno => $seleno_in_contig,
                prot   => $prot
            };
        }
    }

    for (sort keys %visited) {
        my $value = $visited{$_};

        my $status = "not valid";
        $status = $translateStatus{join("", @{$value->{status}})}
          if (exists $value->{status} && defined $value->{status});

        my $exon_string = "$value->{start}-$value->{end}";
        $exon_string = join ",", (map { join "-", @$_ } @{$value->{exons}}) if ($value->{exons});

        print "# alignment\t";
        print join "\t",
          ($exon_string, $status, join(",", @{$value->{ids}}), $value->{prot}, $value->{seleno});
        print "\n";
    }

    my @cds = grep { defined $_->{status} } values %visited;

    # CDS: sort by completeness => sort by length
    my ($current_cds) = sort {
        (completeness_ranking($a->{status}) <=> completeness_ranking($b->{status}))
          || length($b->{prot}) <=> length($a->{prot})

    } @cds;

    unless ($current_cds) {
        return 0;
    } else {
        my $status = $translateStatus{join("", @{$current_cds->{status}})}
          if (defined $current_cds->{status});

        my $exon_string = "$current_cds->{start}-$current_cds->{end}";
        $exon_string = join ",", (map { join "-", @$_ } @{$current_cds->{exons}})
          if ($current_cds->{exons});

        # replace stop codons by one letter code for selenocystein (U)
        $current_cds->{prot} =~ s/\*/U/g if ($current_cds->{seleno});

        my $seleno_string = "NA";
        $seleno_string = join ",",
          (map { $_ * 3 + $current_cds->{start} } @{$current_cds->{seleno}})
          if ($current_cds->{seleno} ne "NA");

        print "alignment\t";
        print join "\t",
          (
            $exon_string, $status, join(",", @{$current_cds->{ids}}),
            $current_cds->{prot}, $seleno_string
          );
        print "\n";
    }

    return 1;
}

sub findIntrons {
    my $seq           = shift;
    my $contig_string = shift;
    my $min_length    = shift || $opt{'intron-length'};

    my @spans;
    my $l = length($seq);
    while ($seq =~ /(-{$min_length,})/g) {
        next if ($-[0] == 0);
        next if ($+[0] == $l);

        my $start_seq = lc(substr($contig_string, $-[0],     2));
        my $stop_seq  = lc(substr($contig_string, $+[0] - 2, 2));

        # conservative: only GT/AG intron start/end
        if ($start_seq eq "gt" && $stop_seq eq "ag") {
            push @spans, [$-[0], $+[0]];
        }
    }

    return \@spans;
}

sub getIntronLess {
    my $string  = shift;
    my @introns = @{shift()};
    my $start   = shift;
    my $stop    = shift;

    # zero based introns, start and stop

    (my $string_seq = $string->seq) =~ s/-//g;

    my @exons;
    my @exon_string;
    my $last_start = $start;
    my $end_of_seq = 0;
    my $length = 0;
    
    for (my $i = 0; $i < @introns; $i++) {
        my $end = $introns[$i]->[0] - 1;
        if ($end > $stop) {
            $end        = $stop;
            $end_of_seq = 1;
        }

        my ($e_start, $e_end) = ($last_start, $end);
        $last_start = $introns[$i]->[1];

        # introns are zero based => +1
        my ($re_start, $re_end) =
          ($string->location_from_column($e_start + 1), $string->location_from_column($e_end + 1));
        $re_start = ($re_start) ? $re_start->start : 1;
        $re_end   = ($re_end)   ? $re_end->end     : 1;

        push @exon_string, substr($string_seq, $re_start - 1, ($re_end - $re_start + 1));
        push @exons, [$re_start, $re_end];
        $length += $re_end - $re_start + 1;

        last if ($end_of_seq);
    }

    unless ($end_of_seq) {
        my ($e_start, $e_end) = ($last_start, $stop);

        # introns are zero based => +1
        my ($re_start, $re_end) =
          ($string->location_from_column($e_start + 1), $string->location_from_column($e_end + 1));
        $re_start = ($re_start) ? $re_start->start   : 1;    
        $re_end   = ($re_end)   ? $re_end->start : 1;    

        push @exon_string, substr($string_seq, $re_start - 1, ($re_end - $re_start + 1));
        push @exons, [$re_start, $re_end];
        $length += $re_end - $re_start + 1;
    }

    # adjust coordinates of last exons in order to include last complete exon
    $exons[$#exons]->[1] -= ($length % 3);

    return (join("", @exon_string), [@exons]);
}

sub getLocation {
    my $start  = shift;
    my $end    = shift;
    my $status = shift;
    my $site   = shift || "start_end";

    my $start_type = "EXACT";
    my $end_type   = "EXACT";

    if ($site =~ /start/ && $status->[0] == 0) {
        $start_type = "BEFORE";
    }
    if ($site =~ /end/ && $status->[1] == 0) {
        $end_type = "AFTER";
    }

    return Bio::Location::Fuzzy->new(
        -start     => $start,
        -end       => $end,
        -start_fuz => $start_type,
        -end_fuz   => $end_type,
        -strand    => 1
    );
}

sub getCDSbyLongestOrf {
    my $command =
      "cat $opt{contig} | getorf -filter -auto -find 2 -noreverse | sizeseq -filter -desc | seqret -filter -first |";

    my $orf = Bio::SeqIO->new(-file => $command)->next_seq;
    unless ($orf) {
        warn("EMBOSS failed.");
        return;
    }

    my $desc = "0 - 0";
    if ($orf->desc =~ /\[(.+?)\]/) {
        $desc = $1;
    }

    my ($start, $end) = split /\s-\s/, $desc;

    print join "\t", "longest_orf", "$start-$end", "partial", "NA", $orf->translate->seq, "NA";

    return 1;
}

sub getCDSbyGenscan {
    my $genscan_file;

    if ($opt{'out-genscan'}) {
        $genscan_file = $opt{'out-genscan'};
    } else {
        (my $fh, $genscan_file) = tempfile(UNLINK => 1);
    }

    unless (-e $genscan_file) {
        my $command = "genscan ".$opt{'genscan-matrix'}." $opt{contig} > $genscan_file 2> /dev/null"; 
        my $failed = system($command);
        if ($failed) {
            warn("Genscan failed.");
        }
    }
    my $genscan = Bio::Tools::Genscan->new(-file => $genscan_file);

    my @predictions;
    eval {
        while (my $gene = $genscan->next_prediction()) {
            push @predictions, $gene;
        }
    };
    if ($@) {
        warn("Could not retrieve genscan result. Skipping genscan!");
        return;
    }

    my @locations;
    for my $gene (@predictions) {
        next
          if ($gene->strand == -1);    # we assume that our contig is orientated correctly
        next unless ($gene->exons);
        my $location = Bio::Location::Split->new();

        my @exons      = $gene->exons;
        my $first_exon = shift @exons;
        if ($first_exon->frame > 0) {
            $first_exon->start($first_exon->start + (3 - $first_exon->frame));
        }
        $location->add_sub_Location($first_exon->location);

        for (@exons) {
            $location->add_sub_Location($_->location);
        }
        push @locations, $location;
    }
    $genscan->close();

    unless (@locations) {
        print "# genscan - no results\n";
        return;
    }

    my $location;
    my $identity = 0;

    # find ORF with highest similirity to ortholog (especially in case of
    # fusion transcript we don't want to use the wrong CDS)
    my $location_index = 0;
    my $orf            = 0;
    if (@locations > 0) {
        for my $loc (@locations) {

            my ($fh,     $file_fasta) = tempfile(UNLINK => 1);
            my ($fh_aln, $file_aln)   = tempfile(UNLINK => 1);

            my $out_seq = $contig_seq->trunc($loc->start, $loc->end);

            my $io = Bio::SeqIO->new(
                -file   => ">" . $file_fasta,
                -format => 'fasta'
            );
            $io->write_seq($out_seq);
            $io->write_seq($compareTo_seq);

            my $o_length = $compareTo_seq->length;
            my $c_length = $out_seq->length;
            my $longest  = $o_length > $c_length ? $o_length : $c_length;

            my $command;
            if ($longest < MAFFT_LENGTH) {
                $command =
                  "mafft --localpair --maxiterate 1000 --quiet $file_fasta > $file_aln 2> /dev/null";
            } else {
                $command = "mafft --quiet $file_fasta > $file_aln 2> /dev/null";
            }
            my $failed = system($command);

            if ($failed) {
                warn("Alignment of genscan predicted CDS and ortholog failed!");
                next;
            }

            my $aln = Bio::AlignIO->new(
                -file   => $file_aln,
                -format => 'fasta'
            )->next_aln;

            close $fh;
            close $fh_aln;

            my $cur_identity = sprintf("%.3f", $aln->overall_percentage_identity('short'));

            if ($cur_identity > $identity) {
                $location_index = $orf;
                $identity       = $cur_identity;
            }

            my $prot = $codon_table->translate($out_seq->seq);
            print join "\t", "# genscan", $loc->start . "-" . $loc->end, $cur_identity, $prot, "NA";
            print "\n";

            $orf++;
        }
    }

    # best matching GENSCAN prediction
    $location = $locations[$location_index];

    my $feature;

    # extract all skipped regions
    my @skipped_bases;
    my @sub_locations  = $location->sub_Location();
    my $first_location = shift @sub_locations;
    my $end            = $first_location->end;
    for (@sub_locations) {
        push @skipped_bases, $contig_seq->trunc($end + 1, $_->start - 1)->seq;
        $end = $_->end;
    }
    my $last_location = pop @sub_locations;

    # skipped bases with canonical splice sites? => intron, otherwise pseudogene
    my $completeness = [0, 0, 0];
    if (@skipped_bases > 0) {
        $completeness->[2] = 1;
        for (@skipped_bases) {
            if ($_ !~ /^GT.+?AG$/) {
                $completeness->[2] = 2;
            }
        }
    }

    # determine completeness (START-codon and STOP-codon)
    # risky to assume completeness, but we know that we are dealing with
    # genscan result
    if (uc($contig_seq->subseq($location->start, $location->start + 2)) eq "ATG") {
        $completeness->[0] = 1;
    }
    if ($codon_table->is_ter_codon($contig_seq->subseq($location->end - 2, $location->end))) {
        $completeness->[1] = 1;
    }

    my $translation;
    if (@skipped_bases > 0) {
        $translation .= $contig_seq->trunc($first_location->start, $first_location->end)->seq;
        for (@sub_locations) {
            $translation .= $contig_seq->trunc($_->start, $_->end)->seq;
        }
        $translation .= $contig_seq->trunc($last_location->start, $last_location->end - 3)->seq;
    } else {
        $translation = $contig_seq->trunc($location->start, $location->end - 3)->seq;
    }

    $translation = "" unless (defined $translation);

    my $location_string = join ",", (map { $_->start . "-" . $_->end } $location->sub_Location);

    my $completeness_string = $translateStatus{join("", @$completeness)};
    print join "\t", "genscan", $location_string,
      $completeness_string, "NA",
      $codon_table->translate($translation), "NA";
    print "\n";

    return 1;
}

# Ranks predict CDS by its completeness
# full length => 5p complete => 3p complete => partial
sub completeness_ranking {
    my ($completeness) = @_;
    if (ref($completeness) ne "ARRAY") {
        return 10;
    }
    if ($completeness->[0] && $completeness->[1]) {    # full length
        return 1;
    } elsif ($completeness->[0] && !$completeness->[1]) {    # 5 prime complete
        return 2;
    } elsif (!$completeness->[0] && $completeness->[1]) {    # 3 prime complete
        return 3;
    }

    return 4;
}

sub max {
    my @length = @_;

    my $max = 0;
    for (@lengths) {
        if ($max < $_) {
            $max = $_;
        }
    }
    return $max;

}

1;

__END__

MAFFT Notes:

    Mafft-linsi not suitable for long sequences (in case of scaffolded TTN it
    terminates with errors [+ takes ages and a lot of memory])
    
    Longest sequences in human refseq:
    
    NM_001267550    109224  
    NM_001256850    104301  
    NM_133378       101520  
    NM_133437       82605   
    NM_133432       82404   
    NM_003319       82029   
    NM_024690       43816   
    NM_182961       27748   
    NM_033071       27439   
    NM_001271223    26925   
    
    => use mafft-linsi for alignments with longest sequences smaller than 30kb
    => auto otherwise

EMBOSS Notes:

    -filter - read from stdin, write to stdout
    -auto - no prompts
    -find 2 - sequence between stops
    -noreverse -  ...
    
GENBANK DEFINITION LINE

    DEFINITION - A concise description of the sequence. Mandatory keyword/one or
    more records.  Example: DEFINITION  Heterocephalus glaber BHMT mRNA
    
