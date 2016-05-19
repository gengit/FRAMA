#!/usr/bin/env perl 

=pcd

=head1 INFO

Martin Bens, bensmartin@gmail.com
2014-07-21

=head1 DESCRIPTION

Ortholog based scaffolding of assembled contigs.

=head1 OPTIONS

=over 8

=item B<-contig> contig.fa

Contig.

=item B<-ortholog> orth.fa

Ortholog.

=item B<-fragments> frag.fa

Fasta file with fragments. If -fragment-list or -blast is not specified, all
entries will be considered for scaffolding. Entries must already be orientated
in that case.

=item B<-blast> blast.txt

Blast result with all possible fragments. Strand will be used to orientate
fragment.

=item B<-fragment-list> list.txt

List of IDs which should be used for scaffolding. Sequences in -fragments must
already be orientated in that case.

=item B<-fragment-overlap> 66

Maximum percentage of overlap between fragments [default: 66]. Pairs exceeding
threshold are further examined for removal based on fragment-identity and
identity to reference.

=item B<-fragment-identity>  98

Minimum identity in overlap between pairs which exceed "-fragment-overlap"
required to keep both fragments [default: 98]. If the identity in overlap is
lower, fragment with higher identity to ortholog is kept.

=item B<-allow-gap-embedded> 200

Mininum required length of gap to consider gapfilling of best hit [default:
200].

=item B<-output> out.fa

(optional) Creates fasta with scaffolded contig.

=item B<-out-combined> combined.aln

(optional) Creates (clustalw) alignment with all fragments.

=item B<-out-final> final.aln

(optional) Creates (clustalw) alignment with remaining fragments after
filtering.

=item B<-cpu> 1

(optional) Use mafft with multiple threads [default: 1].

=item B<-cds> 1 20

(optional) Specify CDS region (1-based). If specified, CDS is added to alignments. [relic
of past].

=back

=cut

use strict;
use warnings;

use FindBin;
use lib "$FindBin::Bin/../lib/perl";

use Getopt::Long;
use Data::Dumper;
use Pod::Usage;
use Bio::DB::Fasta;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::AlignIO::fasta;
use Bio::SimpleAlign;
use Set::IntSpan::Fast;

use File::Path qw(make_path remove_tree);
use File::Spec::Functions;
use File::Temp qw(tempfile tempdir);

use GenbankHelper;
use StringHelper;
use Bio::Range;
use POSIX qw(strftime);

pod2usage( -message => "\n\tNo arguments\n", -verbose => 1 )
  if ( @ARGV == 0 );

use constant MAX_INT => 9**9**9;
use constant MIN_INT => -9**9**9;

use constant QUERY_ID     => 0;
use constant TARGET_ID    => 1;
use constant STRAND       => 7;
use constant TARGET_START => 10;
use constant TARGET_END   => 11;

my $GAP_LENGTH = 10;
my $BBH_EXTEND = 30;

my %opt = (
    'reindex'          => 0,
    'cpu'              => 1,
    'min-contribution' => 20,
    'fragment-overlap' => 5,
    'force-output'     => 0,
    'useId'            => 0
);

GetOptions(
    \%opt,               'contig|c=s',
    'ortholog|r=s',      'cds=i{2}',
    'fragments=s',       'blast=s',
    'fragment-list=s',   'output|o=s',
    'out-combined|oc=s', 'out-final|of=s',
    'cpu=i',             'assembling=s',
    'out-filtered=s',    'fragment-overlap=s',
    'help', 'force-output', 'useId', 
) or pod2usage(1);

pod2usage( -verbose => 2 ) if $opt{help};

$opt{'fragment-overlap'} /= 100;


if ( !$opt{contig} || !$opt{ortholog} || !$opt{fragments} ) {
    die "You must specify -contig, -ortholog and -fragments\n";
}

my $io = Bio::SeqIO->new( -file => $opt{contig}, -format => "fasta" );
die "\nCould not read contig file.\n\n" unless ($io);
my $contig_seq = $io->next_seq;

$io = Bio::SeqIO->new( -file => $opt{ortholog}, -format => "fasta" );
die "\nCould not read ortholog file.\n\n" unless ($io);
my $orth_seq = $io->next_seq;
die "\nNo ortholog sequence found.\n\n" unless ($orth_seq);
my $orth_id = $orth_seq->id;

my $db_fragments =
  Bio::DB::Fasta->new( $opt{fragments}, -reindex => $opt{reindex} );
die "Could not read fragment file" unless ($db_fragments);

my $aln    = Bio::SimpleAlign->new();
my $locseq = Bio::LocatableSeq->new(
    '-seq'        => $orth_seq->seq,
    '-display_id' => $orth_seq->id,
    '-start'      => 1,
    '-end'        => $orth_seq->length,
);
$aln->add_seq($locseq);

#print Dumper $orth_seq;
#print Dumper $locseq;
#print Dumper $aln;

# historical (now only for visualiziation purpose)
my $cds_range;
if ( $opt{cds} ) {
    my ( $start, $end ) = @{ $opt{cds} };

    $cds_range = Bio::Range(
        -start => $start,
        -end   => $end
    );
    $aln->add_seq(
        Bio::LocatableSeq->new(
            -id  => "CDS_" . $orth_seq->id,
            -seq => "-" x ( $start - 1 )
              . $orth_seq->subseq( $start, $end )
              . "-" x ( $orth_seq->length - $end )
        )
    );
}

# alignment output: sequences order => Orth, CDS, Consensus, Fragments
my @order_sequences;
if ( $opt{cds} ) {
    @order_sequences = (
        $orth_seq->id, "CDS_" . $orth_seq->id,
        "Consensus",   "BEST_" . $contig_seq->id
    );
} else {
    @order_sequences =
      ( $orth_seq->id, "Consensus", "BEST_" . $contig_seq->id );
}

# get all fragments
my @fragments;
if ( $opt{blast} ) {    # via blast
    open FH, "<", $opt{blast} or die $!;
    while (<FH>) {
        chomp;
        my @e = split "\t";

        push @fragments, [
            $e[QUERY_ID],    # id
            $e[STRAND] == -1 ? 1 : 0,    # reverse complement?
        ];
    }
    close FH;
} elsif ( $opt{'fragment-list'} ) {    # via ID
    open FH, "<", $opt{blast} or die $!;
    while (<FH>) {
        chomp;
        push @fragments, $_;
    }
    close FH;
} elsif ( $opt{'fragments'} ) {        #
    @fragments = grep { !/^__/ } keys %{ $db_fragments->{offsets} };

    # only thing standing between BioPerl version 1.006901 and 1.006924
    #@fragments = $db_fragments->get_all_primary_ids;
}

# Header
print "###\n";
print "# sequences " . $contig_seq->id . "\n";
print "# fragments " . scalar @fragments . "\n";
print "##\n\n";

unless ( @fragments > 0 ) {
    print "# No fragments found. Nothing to do.\n";
    if ($opt{'force-output'} && $opt{output}) {
        my $io = Bio::SeqIO->new( 
            -file => ">" . $opt{output}, 
            -format => "fasta" 
        );
        $io->write_seq( $contig_seq );
    }
    exit;
}

# using prefix for identification of best hit (BEST) and fragments  (FRAGMENT)
my @fragment_seqs;
for (@fragments) {
    my $fragment_seq;
    my $id;
    if ( defined ref($_) && ref($_) eq "ARRAY" ) {
        $fragment_seq = $db_fragments->get_Seq_by_id( $_->[0] );
        $id = "FRAGMENT_" . $_->[0];
        print "# Not found: $_->[0]\n" unless ($fragment_seq);
        $fragment_seq = $fragment_seq->revcom if ( $fragment_seq && $_->[1] );
    } else {
        $fragment_seq = $db_fragments->get_Seq_by_id($_);
        $id = "FRAGMENT_" . $_;
        print "# Not found: $_\n" unless ($fragment_seq);
    }
    push @fragment_seqs,
      Bio::LocatableSeq->new(
        -id  => "FRAGMENT_" . $_,
        -seq => $fragment_seq->seq
      );
}

unless ( @fragment_seqs == @fragments ) {
    die "\nSequences not found.\n\n";
}
push @fragment_seqs,
  Bio::LocatableSeq->new(
    -id  => "BEST_" . $contig_seq->id,
    -seq => $contig_seq->seq
  );

# assembling overlapping contigs
my ( $fragment_seqs_ref, $cluster2contig );
if ( exists $opt{assembling} ) {
    ( $fragment_seqs_ref, $cluster2contig ) =
      assemble_fragments( \@fragment_seqs );
}

# multiple sequence alignment
if ($fragment_seqs_ref) {
    $aln = align_fragments( $aln, $fragment_seqs_ref );
} else {
    $aln = align_fragments( $aln, \@fragment_seqs );
}

# intermediate output
if ( $opt{'out-combined'} ) {
    my $sorted_aln = sortAlignment( $aln, \@order_sequences );
    my $io = Bio::AlignIO->new(
        -file   => ">" . $opt{'out-combined'},
        -format => "clustalw"
    );
    $io->write_aln($sorted_aln);
}

# set cover
$aln = filterFragments($aln);

# no scaffolding possible
if (   ( $opt{cds} && $aln->num_sequences <= 3 )
    || ( !$opt{cds} && $aln->num_sequences <= 2 ) )
{
    # TODO: write contig IDs of clusters to file
    print "\n# Nothing to do.\n";
    
    if ($opt{'force-output'} && $opt{output}) {
        my $io = Bio::SeqIO->new( 
            -file => ">" . $opt{output}, 
            -format => "fasta" 
        );
        $io->write_seq( $contig_seq );
    }

    exit;
}

# final alignment after filtering
if ( $opt{'out-filtered'} ) {
    my $sorted_aln = sortAlignment( $aln, \@order_sequences );
    my $io = Bio::AlignIO->new(
        -file   => ">" . $opt{'out-filtered'},
        -format => "clustalw"
    );
    $io->write_aln($sorted_aln);
}

my $consensus = prepare_consensus($aln);
my $consensus_seq =
  Bio::LocatableSeq->new( -seq => $consensus, -id => "Consensus" );
$aln->add_seq($consensus_seq);

# remove flanking 'N' and gap character in consensus from alignment
my ( $overhang_left, $overhang_right ) =
  StringHelper::getOverhang( $consensus, "N" );
my $remove_columns = getGapPos($consensus);
if ( $overhang_left > 0 ) {
    push @$remove_columns, [ 0, $overhang_left - 1 ];
}
if ( $overhang_right > 0 ) {
    push @$remove_columns,
      [ length($consensus) - $overhang_right, length($consensus) ];
}

if ( @$remove_columns > 0 ) {
    $aln = $aln->remove_columns(@$remove_columns);
}

if ( $opt{output} ) {
    my $io = Bio::SeqIO->new( -file => ">" . $opt{output}, -format => "fasta" );
    my $seq = $aln->get_seq_by_id("Consensus");
    if ($opt{'useId'}) {
        $seq->id($contig_seq->id);
    }
    $io->write_seq( $seq );
}

if ( $opt{'out-final'} ) {
    my $sorted_aln = sortAlignment( $aln, \@order_sequences );
    my $io = Bio::AlignIO->new(
        -file   => ">" . $opt{'out-final'},
        -format => "clustalw"
    );
    $io->write_aln($sorted_aln);
}

# get fragment positions
my %fragment_pos;
my $consensus_string;
for ( $aln->each_seq ) {
    $consensus_string = $_->seq if ( $_->id eq "Consensus" );

    my ( $start, $end ) = StringHelper::getRange( $_->seq, "-" );
    $fragment_pos{ $_->id } = [ $start, $end ];
}

my %map = map { $_ => 1 } @order_sequences;
my @remaining = grep { !$map{$_} } map { $_->id } $aln->each_seq;
push @order_sequences, @remaining;

@fragments = grep { exists $fragment_pos{$_} } @remaining;
@fragments =
  sort { $fragment_pos{$a}->[0] <=> $fragment_pos{$b}->[0] } @fragments;

print join "\t", "best", $contig_seq->id,
  $fragment_pos{ "BEST_" . $contig_seq->id }->[0],
  $fragment_pos{ "BEST_" . $contig_seq->id }->[1];
print "\n";

for (@fragments) {
    ( my $clean_id = $_ ) =~ s/FRAGMENT_//;
    print join "\t", "fragment", $clean_id, $fragment_pos{$_}->[0],
      $fragment_pos{$_}->[1];
    print "\n";
}

while ( $consensus_string =~ /(N{5,})/ig ) {
    print join "\t", "gap", "", $-[0] + 1, $+[0];
    print "\n";
}

if ( $opt{'out-contig'} ) {
    my $io =
      Bio::SeqIO->new( -file => ">" . $opt{'out-contig'}, -format => "fasta" );
    $io->write_seq(
        Bio::Seq->new( -id => $contig_seq->id, -seq => $consensus_string ) );
}

sub count_seq {
    my ( $aln, $pos ) = @_;

    my @count;
    my $counter = 0;
    for my $seq ( $aln->each_seq ) {
        my $letter = substr( $seq->seq, $pos, 1 );

        return [$counter] if ( $seq->id =~ /^BEST/ && $letter ne "X" );

        push @count, $counter if ( $letter ne "X" );
        $counter++;
    }

    if ( @count > 0 ) {
        return \@count;
    } else {
        return [-1];
    }
}

sub prepare_consensus {
    my ($aln) = @_;

    # prepare alignment; replace overhanging "-" with "X"
    my $consensus_aln = Bio::SimpleAlign->new();
    for ( $aln->each_seq ) {
        next if ( $_->id =~ /$orth_id/ );

        my $seq_string = StringHelper::replaceOverhangChar( $_->seq, "-", "X" );
        my $locseq = Bio::LocatableSeq->new(
            -seq => $seq_string,
            -id  => $_->id
        );
        $consensus_aln->add_seq($locseq);
    }

    # get ranges and identities of all contigs
    my @ranges;
    my @ids;
    my $length = $consensus_aln->length_flushed;
    for ( $consensus_aln->each_seq ) {
        my ( $start, $end ) = StringHelper::getRange( $_->seq, "X" );
        my $range = Bio::Range->new( -start => $start, -end => $end );
        push @ranges, $range;
        push @ids,    $_->id;
    }

    # get number of sequences for each position (to identify overlapping regions)
    my @num_sequences = map { 0 } 0 .. ( $length - 1 );
    for ( 0 .. ( $length - 1 ) ) {
        $num_sequences[$_] = count_seq( $consensus_aln, $_ );
    }

    my @e = grep { defined ref($_) and ref($_) eq "ARRAY" and @{$_} > 2 }
      @num_sequences;
    if ( @e > 0 ) {
        print STDERR
          "WARNING: More than 2 overlapping sequences not implemented yet.\n";
    }

    # get overlapping regions as [Bio::Range, [pos sequence1, pos sequence2]]
    my $start = 0;
    my $last  = $num_sequences[0];
    if ( @$last > 1 ) {
        $start = 1;
    }
    $last = join( " ", @$last );

    my @overlapping;
    for ( my $i = 1; $i < @num_sequences; $i++ ) {
        my $index = $num_sequences[$i];
        my $current = join( " ", @$index );

        if ( @$index > 1 ) {
            if ( $current ne $last ) {
                if ( $start > 0 ) {
                    push @overlapping,
                      [
                        Bio::Range->new( -start => $start, -end => $i ),
                        [ split " ", $last ]
                      ];
                }
                $start = $i + 1;
            }
        } else {
            if ( $start > 0 ) {
                push @overlapping,
                  [
                    Bio::Range->new( -start => $start, -end => $i ),
                    [ split " ", $last ]
                  ];
            }
            $start = 0;
        }
        $last = $current;
    }

    if ($start) {
        push @overlapping,
          [
            Bio::Range->new( -start => $start, -end => scalar @num_sequences ),
            [ split " ", $last ]
          ];
    }

    for my $overlap (@overlapping) {

        if ( $overlap->[0]->length > 5 ) {
            my $first_half_end =
              $overlap->[0]->start + int( $overlap->[0]->length / 2 ) - 1;

            my $first_range = Bio::Range->new(
                -start => $overlap->[0]->start,
                -end   => $first_half_end
            );
            my $second_range = Bio::Range->new(
                -start => $first_half_end + 1,
                -end   => $overlap->[0]->end
            );

            # first half
            my $first_ID1 =
              getIdentitySlice( $aln, $orth_id, $ids[ $overlap->[1]->[0] ],
                $first_range, 0 );
            my $first_ID2 =
              getIdentitySlice( $aln, $orth_id, $ids[ $overlap->[1]->[1] ],
                $first_range, 0 );
            my $first_index = $first_ID1 > $first_ID2 ? 0 : 1;

            my $second_ID1 =
              getIdentitySlice( $aln, $orth_id, $ids[ $overlap->[1]->[0] ],
                $second_range, 0 );
            my $second_ID2 =
              getIdentitySlice( $aln, $orth_id, $ids[ $overlap->[1]->[1] ],
                $second_range, 0 );
            my $second_index = $second_ID1 > $second_ID2 ? 0 : 1;

            for ( $first_range->start - 1 .. $first_range->end - 1 ) {
                $num_sequences[$_] = [ $overlap->[1]->[$first_index] ];
            }
            for ( $second_range->start - 1 .. $second_range->end - 1 ) {
                $num_sequences[$_] = [ $overlap->[1]->[$second_index] ];
            }

            print "##\n";
            print "# Overlap between "
              . $ids[ $overlap->[1]->[0] ] . " and "
              . $ids[ $overlap->[1]->[1] ] . "\n";
            print "#    Range: "
              . $overlap->[0]->start . "-"
              . $overlap->[0]->end . "\n";
            print "#    Identity 1st half: $first_ID1 - $first_ID2; "
              . $first_range->start . "-"
              . $first_range->end . " => "
              . $ids[ $overlap->[1]->[$first_index] ] . "\n";
            print "#    Identity 2nd half: $second_ID1 - $second_ID2; "
              . $second_range->start . "-"
              . $second_range->end . " => "
              . $ids[ $overlap->[1]->[$second_index] ] . "\n";
            print "\n";
        } else {
            my $first_ID1 =
              getIdentitySlice( $aln, $orth_id, $ids[ $overlap->[1]->[0] ],
                $overlap->[0], 0 );
            my $first_ID2 =
              getIdentitySlice( $aln, $orth_id, $ids[ $overlap->[1]->[1] ],
                $overlap->[0], 0 );
            my $first_index = $first_ID1 > $first_ID2 ? 0 : 1;
            for ( $overlap->[0]->start - 1 .. $overlap->[0]->end - 1 ) {
                $num_sequences[$_] = [ $overlap->[1]->[$first_index] ];
            }
            print "##\n";
            print "# Overlap between "
              . $ids[ $overlap->[1]->[0] ] . " and "
              . $ids[ $overlap->[1]->[1] ] . "\n";
            print "#    Range: "
              . $overlap->[0]->start . "-"
              . $overlap->[0]->end . "\n";
            print "#    Identity: $first_ID1 - $first_ID2; "
              . $overlap->[0]->start . "-"
              . $overlap->[0]->end . " => "
              . $ids[ $overlap->[1]->[$first_index] ] . "\n";
            print "\n";
        }
    }

    my @chars;
    my $counter = 0;
    for (@num_sequences) {
        my $pos = $_->[0];
        if ( $pos == -1 ) {
            push @chars, "N";
        } else {
            push @chars,
              substr( $consensus_aln->get_seq_by_pos( $pos + 1 )->seq,
                $counter, 1 );
        }
        $counter++;
    }

    return join "", @chars;
}

sub set_cover {
    my ( $S, $R, $w ) = @_;
    my $minCost    = MAX_INT;
    my $minElement = -1;
    for ( my $i = 0; $i < @$S; $i++ ) {
        my $i_set = $S->[$i]->intersection($R);

        next unless ( $i_set->cardinality > $opt{'min-contribution'} );

        # TODO: sqrt
        my $cost = $w->[$i] / $i_set->cardinality;

        if ( $cost < $minCost ) {
            $minCost    = $cost;
            $minElement = $i;
        }
    }
    return undef if ( $minElement == -1 );

    return ( $minElement, $minCost );
}

sub filterFragments {
    my ($aln) = @_;

    my %output;
    my ( @ranges, @weights, @ids );
    for my $seq ( $aln->each_seq ) {
        ( my $clean_id = $seq->id ) =~ s/(.+?)\/(.+)/$1/g;
        next unless ( $clean_id =~ /^FRAGMENT_/ || $clean_id =~ /^BEST/ );

        push @ids, $clean_id;

        # range
        my ( $start, $end ) = StringHelper::getRange( $seq->seq, "-" );
        my $range = Set::IntSpan::Fast->new("$start-$end");
        push @ranges, $range;

        my $identity =
          getIdentitySlice( $aln, $orth_id, $seq->id, [ $start, $end ], 1 );
        my $weight =
          ( $clean_id =~ /^BEST/ ) ? 0 : sprintf( "%.0f", ( 100 - $identity ) );
        push @weights, $weight;

        # ouptut [identity, length]
        $output{$clean_id} = [ $identity, $range->cardinality, $weight ];
    }

    # output
    print "##\n";
    print "# Fragments\n";
    for ( sort keys %output ) {
        print "#    $_ "
          . $output{$_}->[0] . "; "
          . $output{$_}->[1] . "; "
          . $output{$_}->[2] . "\n";
    }
    print "\n";

    my $length         = $aln->length_flushed;
    my $alignment_span = Set::IntSpan::Fast->new("1-$length");

    # flag resulting fragments
    my @result = map { 0 } @ranges;

    my ( @order_frag, @cost );
    while () {
        my ( $S_i, $cost ) = set_cover( \@ranges, $alignment_span, \@weights );

        # no contig found
        last if ( not defined $S_i );

        $result[$S_i] = 1;
        push @order_frag, $S_i;
        push @cost,       $cost;

        # remove fragment from
        $alignment_span->remove_from_string( $ranges[$S_i]->as_string );
    }

    my $count = 0;
    my ( @ranges_bio, %ids, @ids_bio );
    for my $seq ( $aln->each_seq ) {
        if ( $seq->id =~ /^FRAGMENT_/ || $seq->id =~ /^BEST/ ) {
            if ( $result[$count] ) {
                my ( $start, $end ) = StringHelper::getRange( $seq->seq, "-" );
                push @ranges_bio,
                  Bio::Range->new( -start => $start, -end => $end );
                $ids{ $seq->id } = 1;
                push @ids_bio, $seq->id;
            }
            $count++;
        }
    }

    # remove shorter contigs embedded in longer contigs
    my @range_order = sort {
             $ranges_bio[$a]->start <=> $ranges_bio[$b]->start
          || $ranges_bio[$b]->length <=> $ranges_bio[$a]->length
    } 0 .. $#ranges_bio;

    my %completely_embedded;
    my @indices;
    for ( my $i = 0; $i < @ranges_bio; $i++ ) {
        my $id_A = $ids_bio[ $range_order[$i] ];
        next if ( $id_A =~ /^BEST/ );
        next if ( $completely_embedded{$id_A} );

        my $range_A = $ranges_bio[ $range_order[$i] ];

        for ( my $j = $i + 1; $j < @ranges_bio; $j++ ) {
            my $id_B = $ids_bio[ $range_order[$j] ];
            next if ( $id_B =~ /^BEST/ );
            next if ( $completely_embedded{$id_A} );

            my $range_B = $ranges_bio[ $range_order[$j] ];

            if ( $range_A->contains($range_B) ) {
                $completely_embedded{$id_B} = 1;
                $ids{$id_B}                 = 0;
                push @indices, $range_order[$j];
            }
        }
    }

    my $filtered_aln = Bio::SimpleAlign->new();
    $filtered_aln->add_seq( $aln->get_seq_by_id($orth_id) );

    print "##\n";
    print "# Keeping sequences\n";
    for ( $aln->each_seq ) {
        if ( $ids{ $_->id } ) {
            $filtered_aln->add_seq($_);
            print "#    " . $_->id . "\n";
        }
    }
    print "\n";

    die "Filtering went wrong. This should not happend\n"
      unless ($filtered_aln);
    $filtered_aln = $filtered_aln->remove_columns( ['all_gaps_columns'] );

    return $filtered_aln;
}

sub overhangTrim {
    my ( $aln, $refID, $refPos ) = @_;

    my $ref_seq;
    if ($refID) {
        $ref_seq = $aln->get_seq_by_id($refID);
    } elsif ( $refPos > 0 ) {
        $ref_seq = $aln->get_seq_by_pos($refPos);
    } else {
        $ref_seq = $aln->get_seq_by_pos(1);
        print STDERR
          "Not reference sequence specified. Using first sequence in alignment: "
          . $ref_seq->id . "\n";
    }

    my $start = 1;
    my $end   = $aln->length_flushed;

    my $left_overhang  = 0;
    my $right_overhang = 0;

    my $new_aln = $aln;

    if ( $ref_seq->seq =~ /^(-+)/ ) {
        $left_overhang = length($1);
        $start         = length($1) + 1;
    }
    if ( $ref_seq->seq =~ /(-+)$/ ) {
        $right_overhang = length($1);
        $end -= length($1);
    }

    # get slice of aligned part
    if ( $left_overhang != 0 || $right_overhang != 0 ) {
        $new_aln = $new_aln->slice( $start, $end )
          if ( $left_overhang > 1 && $right_overhang > 1 );
    }

    return $new_aln;
}

# Computes pairwise identity between to specified sequences (IDs) within a
# specific region (Bio::Range). Sequences must be aligned and provided in
# as Bio::SimpleAlign.
#
# Gaps longer than 9bp are removed beforehand.
#
sub getIdentitySlice {
    my $aln      = shift;
    my $id1      = shift;
    my $id2      = shift;
    my $range    = shift;
    my $clipping = shift || 0;

    if ( defined ref($range) && ref($range) eq "ARRAY" ) {
        $range = Bio::Range->new( -start => $range->[0], -end => $range->[1] );
    }
    my $seq1 = $aln->get_seq_by_id($id1)->trunc( $range->start, $range->end );
    my $seq2 = $aln->get_seq_by_id($id2)->trunc( $range->start, $range->end );

    my $new_aln = Bio::SimpleAlign->new();
    $new_aln->add_seq($seq1);
    $new_aln->add_seq($seq2);
    $new_aln = $new_aln->remove_columns( ['all_gaps_columns'] );

    if ($clipping) {
        $new_aln = overhangTrim( $new_aln, $id1 );
        $new_aln = overhangTrim( $new_aln, $id2 );
    }

    return sprintf "%.3f", $new_aln->overall_percentage_identity('align');
}

#
#  Returns consensus string. Consensus character has to be present in at least
#  50% of all sequences or "N" is used as consensus character. Ties are broken
#  by alphabetical order of nucleotides. Character 'X' does not contribute to
#  the number of sequences, meaning that for a position with two characters "A"
#  and "X", "A" is used as consensus character. 'X'-only columns are represent
#  by 'N'. If number of gaps equals number of (different) letters, "N" will be
#  used. If number of gaps equals consensus character, consensus character will
#  be printed in lowercase.
#
#  Consensus is compute between fragments only. Sequence of best hit is used
#  completely.
#

sub consensus {
    my $aln = shift;

    my $out = "";
    for ( 0 .. ( $aln->length_flushed - 1 ) ) {
        $out .= consensus_position( $aln, $_ );
    }
    return $out;
}

sub consensus_position {
    my $aln = shift;
    my $pos = shift;

    my $threshold_pct = 50;

    my %letters;
    $letters{"X"} = 0;
    for my $seq ( $aln->each_seq ) {
        my $letter = substr( $seq->seq, $pos, 1 );

        return $letter if ( $seq->id =~ /^BEST/ && $letter ne "X" );

        $letters{$letter}++;
    }

    my $letter = 'N';

    my $number_of_sequences = $aln->num_sequences();

    # fill space between non overlapping fragments with N
    if ( $number_of_sequences == $letters{"X"} ) {
        return $letter;
    }
    my $threshold =
      ( $number_of_sequences - $letters{"X"} ) * $threshold_pct / 100.;
    my $count = -1;

    my @possible_letters;    # collect letters of same high score
    foreach my $key ( sort keys %letters ) {
        next if ( $key eq "X" );
        my $num = $letters{$key};

        if ( $letters{$key} >= $threshold ) {
            if ( $letters{$key} > $count ) {
                $letter           = $key;
                $count            = $letters{$key};
                @possible_letters = ($letter);
            } elsif ( $letters{$key} == $count ) {
                push @possible_letters, $key;
            }
        }
    }

    # first non gap character in alphabetal order
    if ( @possible_letters > 1 ) {
        ($letter) = grep { $_ ne "-" } sort @possible_letters;
    }

    return $letter;
}

sub getGapPos {
    my $seq = shift;
    my $one_based = shift || 0;

    my @spans;
    while ( $seq =~ /(-+)/g ) {
        push @spans, [ $-[0], $+[0] ];
    }

    if ($one_based) {
        for (@spans) {
            $_->[0]++;
        }
    } else {
        for (@spans) {
            $_->[1]--;
        }
    }

    return \@spans;
}

sub getCoordinatesCharacterSpan {
    my $seq       = shift;
    my $char      = shift;
    my $span      = shift;
    my $one_based = shift || 0;

    my $regexp = "(" . $char . "{$span,}" . ")";

    my @spans;
    while ( $seq =~ /$regexp/g ) {
        push @spans, [ $-[0], $+[0] ];
    }

    if ($one_based) {
        for (@spans) {
            $_->[0]++;
        }
    } else {
        for (@spans) {
            $_->[1]--;
        }
    }

    return @spans;

}

sub align_fragments {
    my $aln       = shift;
    my $fragments = shift;

    my ( $fh_a, $msa_file ) = tempfile();
    my $align_io =
      Bio::AlignIO->new( -file => ">" . $msa_file, -format => 'fasta' );
    $align_io->write_aln($aln);

    my ( $fh_b, $fragments_file ) = tempfile();
    my $frag_io =
      Bio::SeqIO->new( -file => ">" . $fragments_file, -format => 'fasta' );
    for (@$fragments) {
        $frag_io->write_seq($_);
    }

    my ( $fh_o, $output_file ) = tempfile();

    my $command =
      "mafft --thread $opt{cpu} --addfragments $fragments_file $msa_file > $output_file 2> /dev/null  ";
    system($command);

    my $combined =
      Bio::AlignIO->new( -file => $output_file, -format => 'fasta' );

    #my $combined = Bio::SeqIO->new(-file => $output_file, -format => 'fasta');

    unless ($combined) {
        print "# Error executing:\n";
        print "# \t$command\n";

        return undef;

        close($fh_a);
        close($fh_b);
    }
    my $alignment = $combined->next_aln;

    #print Dumper $alignment;

    unless ($alignment) {
        die
          "Sequences for scaffolding could not be aligned. Check path and workability of mafft.\n";
    }

    return $alignment;
}

sub sortAlignment {
    my $aln   = shift;
    my $order = shift;

    my $num = 0;
    my %order;
    for (@$order) {
        $order{$_} = ++$num;
    }
    my ( @seqs, @ids );
    for my $seq ( $aln->each_seq ) {
        push @seqs, $seq;
        push @ids,  $seq->display_id;
        unless ( exists $order{ $seq->display_id } ) {
            $order{ $seq->display_id } = ++$num;
        }
    }
    my @sorted_seq =
      map { $_->[1] }
      sort { $a->[0] <=> $b->[0] } map { [ $order{ $_->id() }, $_ ] } @seqs;
    my $newaln = Bio::SimpleAlign->new();
    for (@sorted_seq) {
        $newaln->add_seq($_);
    }

    return $newaln;
}

sub Bio::SimpleAlign::length_flushed {
    my $self = shift;
    return $self->get_seq_by_pos(1)->length;
}

sub Bio::AlignIO::fasta::next_aln {
    my $self = shift;
    my ($width) = $self->_rearrange( [qw(WIDTH)], @_ );
    $self->width( $width || 60 );

    my (
        $start, $end,      $name,     $seqname, $seq,  $seqchar,
        $entry, $tempname, $tempdesc, %align,   $desc, $maxlen
    );
    my $aln = Bio::SimpleAlign->new();

    while ( defined( $entry = $self->_readline ) ) {
        chomp $entry;
        if ( $entry =~ s/^>\s*(\S+)\s*// ) {
            $tempname = $1;
            chomp($entry);
            $tempdesc = $entry;
            if ( defined $name ) {
                $seqchar =~ s/\s//g;
                if ( $name =~ /(\S+)\/(\d+)-(\d+)/ ) {
                    $seqname = $1;
                    $start   = $2;
                    $end     = $3;
                } else {
                    $seqname = $name;
                    $start   = 1;
                    $end     = $self->_get_len($seqchar);
                }

                $seq = Bio::LocatableSeq->new(
                    '-seq'         => $seqchar,
                    '-display_id'  => $seqname,
                    '-description' => $desc,
                    '-start'       => $start,
                    '-end'         => $end,
                    '-alphabet'    => $self->alphabet,
                );
                $aln->add_seq($seq);
                $self->debug("Reading $seqname\n");
            }
            $desc    = $tempdesc;
            $name    = $tempname;
            $desc    = $entry;
            $seqchar = "";
            next;
        }

        # removed redundant symbol validation
        # this is already done in Bio::PrimarySeq
        $seqchar .= $entry;
    }

    #  Next two lines are to silence warnings that
    #  otherwise occur at EOF when using <$fh>
    $name    = "" if ( !defined $name );
    $seqchar = "" if ( !defined $seqchar );
    $seqchar =~ s/\s//g;

    #  Put away last name and sequence
    if ( $name =~ /(\S+)\/(\d+)-(\d+)$/ ) {
        $seqname = $1;
        $start   = $2;
        $end     = $3;
    } else {
        $seqname = $name;
        $start   = 1;
        $end     = $self->_get_len($seqchar);
    }

    # This logic now also reads empty lines at the
    # end of the file. Skip this is seqchar and seqname is null
    unless ( length($seqchar) == 0 && length($seqname) == 0 ) {
        $seq = Bio::LocatableSeq->new(
            -seq         => $seqchar,
            -display_id  => $seqname,
            -description => $desc,
            -start       => $start,
            -end         => $end,
            -alphabet    => $self->alphabet,
        );
        $aln->add_seq($seq);
        $self->debug("Reading $seqname\n");
    }
    my $alnlen = $aln->length;
    foreach my $seq ( $aln->each_seq ) {
        if ( $seq->length < $alnlen ) {
            my ($diff) = ( $alnlen - $seq->length );
            $seq->seq( $seq->seq() . "-" x $diff );
        }
    }

    # no sequences means empty alignment (possible EOF)
    return $aln if $aln->num_sequences;
    return;
}

1;

__END__

