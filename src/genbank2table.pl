#!/usr/bin/env perl

=pod

=head1 INFO

Martin Bens, bensmartin@gmail.com
12/16/2013 17:01:51

=head1 DESCRIPTION

Create one line for each GenBank entry.

=head1 OPTIONS

=over 8

=item B<-assembly>

Genbank file created with pipeline. "N" characters are removed from total
length of transcript and CDS.

=item B<-genbank> transcript.gb

Input file.

=item B<-output> file.csv

Output file. [default: STDOUT]

=item B<-cds>

CDS coordinates based on mRNA sequence (ID, Sym, start, end)

=back

=cut

use strict;
use warnings;

use Getopt::Long;
use File::Basename;
use Data::Dumper;
use Pod::Usage;
use Bio::SeqIO;

use constant CDS_TOOLS => qw(genscan alignment longest_orf);

pod2usage( -message => "\n\tNo arguments\n", -verbose => 1 ) if ( @ARGV == 0 );

my $man  = 0;
my $help = 0;
my ( $genbank_file, $output );
my $assembly = 0;
my $cds      = 0;
my $mrna     = 0;
my $tag      = "mRNA";
GetOptions(
    'assembly'  => \$assembly,
    'output=s'  => \$output,
    'cds'       => \$cds,
    'genbank=s' => \$genbank_file,
    'help|?'    => \$help,
    'mrna'      => \$mrna,
    'man|m'     => \$man,
    'tag=s'     => \$tag,
) or pod2usage(1);

pod2usage(1) if $help;
pod2usage( -verbose => 2 ) if $man;

my $io = Bio::SeqIO->new( -file => $genbank_file, -format => 'genbank' );
my $builder = $io->sequence_builder;
$builder->want_none();
$builder->add_wanted_slot( 'features', 'seq' );  # we need sequence to determine
                                                 # number of 'N' characters

my $out;
if ($output) {
    open $out, ">", $output or die $!;
} else {
    $out = *STDOUT;
}

if ($assembly) {
    assembly();
} elsif ($cds) {
    cdsOnly();
} else {
    regular_genbank();
}

sub cdsOnly {
    print $out join "\t", "#accession", "symbol", "start", "end";
    print $out "\n";

    while ( my $seq = $io->next_seq ) {
        next if ( $seq->molecule ne "mRNA" );
        my ($gene) = $seq->get_SeqFeatures('mRNA');
        my $sym;
        if ($gene) {
            ($sym) = $gene->get_tag_values('gene')
              if ( $gene && $gene->has_tag('gene') );
            ($sym) = $gene->get_tag_values('note')
              if ( $gene && !defined $sym && $gene->has_tag('note') );
            $sym = "NA" unless ($sym);
        } else {
            my ($gene) = $seq->get_SeqFeatures('CDS');
            ($sym) = $gene->get_tag_values('gene')
              if ( $gene && $gene->has_tag('gene') );
            ($sym) = $gene->get_tag_values('note')
              if ( $gene && !defined $sym && $gene->has_tag('note') );
            ($sym) = $gene->get_tag_values('product')
              if ( $gene && !defined $sym && $gene->has_tag('product') );
            $sym = "NA" unless ($sym);
        }

        my ($cds) = $seq->get_SeqFeatures('CDS');

        if ($cds) {
            my $start = $cds->start;
            if ( $gene && $gene->start > 1 ) {
                $start = $cds->start - $gene->start + 1;
            }

            my $end = $start + $cds->spliced_seq->length - 1;

            print $out join "\t", $seq->id, $sym, $start, $end;
            print $out "\n";
        } else {
            print STDERR $seq->id . " no CDS found\n";
        }
    }
    close $out;
}

sub assembly {
    print $out join(
        "\t",
        (
            "assembly",      "hit_type",
            "ortholog",      "symbol",
            "length",        "cds_length",
            "num_fragments", "fragments",
            "cds_status",    "cds_prediction_tool",
            "clipped",       "contamination"
        )
    ) . "\n";
    while ( my $seq = $io->next_seq ) {

        # sequence length (excluding N)
        my $seq_string        = $seq->seq;
        my $count             = ( $seq_string =~ tr/N// );
        my $transcript_length = length($seq_string) - $count;

        # source: hit type, ortholog
        my ($source) = $seq->get_SeqFeatures('source');
        unless ($source) {
            print STDERR $seq->id . " has no 'source' feature. Skipping.\n";
            next;
        }
        my ($note) = $source->get_tag_values('note');
        my ( $hit_type, undef, $orth_id ) = ( "", undef, "" );
        if ($note) {
            ( $hit_type, undef, $orth_id ) = split( ":", $note );
        }

        # mRNA: symbol, clipping
        my ($gene) = $seq->get_SeqFeatures('mRNA');
        unless ($gene) {
            print STDERR $seq->id . " has no 'mRNA' feature. Skipping.\n";
            next;
        }

        # number of clipped bases
        my ( undef, $clipped ) = split ":", ( $gene->get_tag_values('note') )[0]
          if ( $gene->has_tag('note') );
        if ($clipped) {
            $clipped = $source->end - $gene->end;
        }

        ($gene) = $gene->get_tag_values('gene');

        # CDS (TODO: should check for molecule first)
        my ( $cds_length, $cds_note, $cds_prediction_tool ) = ( "NA", "NA", "NA" );
        my ($cds) = $seq->get_SeqFeatures('CDS');
        unless ($cds) {
            print STDERR $seq->id . " has no 'CDS' feature. Skipping.\n";
            next;
        }
        ($cds_note) = $cds->get_tag_values('note') if ( $cds->has_tag('note') );

        if ( $cds->has_tag('inference') ) {
            ($cds_prediction_tool) = $cds->get_tag_values('inference');
            for (CDS_TOOLS) {
                if ( $cds_prediction_tool =~ /$_/ ) {
                    $cds_prediction_tool = $_;
                    last;
                }
            }
        }

        # CDS length (excluding N)
        $seq_string = $cds->spliced_seq->seq;
        $count      = ( $seq_string =~ tr/N// );
        $cds_length = length($seq_string) - $count;

        # misc_features (scaffolding, contamination)
        my @fragments = $seq->get_SeqFeatures('misc_feature');
        my @frag_ids;
        my $contamination = "NA";
        for (@fragments) {
            next unless ( $_->has_tag('note') );
            my ($note) = $_->get_tag_values('note');

            if (   $note =~ /^best contig/
                || $note =~ /sequence inferred during scaffolding/ )
            {
                ( undef, my $contig ) = split( ":", $note );
                push @frag_ids, $contig if ($contig);
            } elsif ( $note =~ /(.+?):(.+?):(.+?)/ ) {

                # contamination
                if ( $1 eq "terminal" || $1 eq "internal" ) {
                    $contamination = $note;
                }
            }
        }

        my $num_frag = @frag_ids > 0 ? scalar @frag_ids : 0;
        my $frag_ids = @frag_ids > 0 ? join( ",", @frag_ids ) : "NA";

        print $out join(
            "\t",
            (
                $seq->id,           $hit_type,
                $orth_id,           $gene,
                $transcript_length, $cds_length,
                $num_frag,          $frag_ids,
                $cds_note,          $cds_prediction_tool,
                $clipped,           $contamination
            )
        );
        print $out "\n";
    }
    close $out;
}

sub regular_genbank {
    print $out join(
        "\t",
        (
            "accession",    "nucl_version", "length", "cds_length",
            "symbol",       "geneid",       "ccds",   "protein",
            "prot_version", "organism"
        )
    ) . "\n";
    while ( my $seq = $io->next_seq ) {
        my ($gene) = $seq->get_SeqFeatures('CDS');
        next unless ($gene);

        my $sym = "NA";
        if ( $gene->has_tag('gene') ) {
            ($sym) = $gene->get_tag_values('gene');
        } elsif ( $gene->has_tag('note') ) {

            # sometimes symbol seem to be saved as note..
            ($sym) = $gene->get_tag_values('note');
        } elsif ( $gene->has_tag('product') ) {
            ($sym) = $gene->get_tag_values('product');
        }

        my $proteinid = "NA";
        my $version   = "NA";
        if ( $gene->has_tag('protein_id') ) {
            ($proteinid) = $gene->get_tag_values('protein_id');
            ( $proteinid, $version ) = split "\\.", $proteinid
              if ( $proteinid && $proteinid =~ /\./ );
        }

        my $geneid;
        my $ccds = "NA";
        if ( $gene->has_tag('db_xref') ) {
            ($geneid) =
              map { /GeneID:(.+)/; $1 }
              grep { /GeneID:(.+)/ } $gene->get_tag_values('db_xref');
            ($geneid) =
              map { /GI:(.+)/; $1 }
              grep { /GI:(.+)/ } $gene->get_tag_values('db_xref')
              unless ($geneid);
            ($ccds) =
              map { /CCDS:(.+?)\./; $1 }
              grep { /CCDS:(.+)/ } $gene->get_tag_values('db_xref');
        }

        $proteinid = "NA" unless ($proteinid);
        $version   = "NA" unless ($version);
        $ccds      = "NA" unless ($ccds);
        $geneid    = "NA" unless ($geneid);

        my ($source) = $seq->get_SeqFeatures('source');
        my ($org)    = $source->has_tag('organism') ? $source->get_tag_values('organism') : ( "NA" );

        my $seq_version = $seq->version ? $seq->version : "NA";

        print $out join(
            "\t",
            (
                $seq->id, $seq_version, $seq->length,
                $gene->end - $gene->start + 1,
                $sym, $geneid, $ccds, $proteinid, $version, $org
            )
        ) . "\n";
    }
    close $out;
}

1;

__END__


