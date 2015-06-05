#!/usr/local/bin/perl

use strict;
use warnings;

package GenbankHelper;

use Bio::Index::GenBank;

sub getSymbol {
    my $seq = shift;

    my $gene_symbol;
    eval {

        my ($tag) = grep { $_->primary_tag eq 'CDS' } $seq->get_SeqFeatures;
        ($tag) = grep { $_->primary_tag eq 'gene' } $seq->get_SeqFeatures unless (defined $tag);
        return undef unless ( defined $tag );

        # looking for gene symbol in gene -> note -> product;
        # takes any text
        if ( $tag->has_tag('gene') ) {
            ($gene_symbol) = $tag->get_tag_values('gene');
        } elsif ($tag->has_tag('note')) {
            ($gene_symbol) = $tag->get_tag_values('note');
        } elsif ($tag->has_tag('product')) {
            ($gene_symbol) = $tag->get_tag_values('product');
        } else {
            warn "No symbol assigned!\n";
            return undef;
        }
    };

    return $gene_symbol;
}

sub get_dbxref {
    my $seq  = shift;
    my $type = shift;

    $type = "GeneID" unless ( defined $type );

    my ($tag) = grep { $_->primary_tag eq 'gene' } $seq->get_SeqFeatures;
    return undef unless ( defined $tag );

    my $gene_symbol;
    if ( $tag->has_tag('db_xref') ) {
        my @db_xrefs = $tag->get_tag_values('db_xref');
        for (@db_xrefs) {
            if (/$type:(.+)/) {
                $gene_symbol = $1;
            }
        }
    }

    return $gene_symbol;
}

sub getIndex {
    my $file = shift;

    unless ($file) {
        warn("Genbank file not defined!\n");
        return undef;
    }
    unless (-e $file) {
        warn("Genbank file not found: $file\n");
        return undef;
    }


    my $index_file = $file . ".indexdb";
    my $db_orth;

    if (@_) {
        my %hash = @_;
        if ($hash{'-reindex'} == 1 && -e $index_file) {
            unlink($index_file);
        }
    }

    eval {
        #$db_orth = Bio::Index::GenBank->new();
        $db_orth = Bio::Index::GenBank->new( -filename => $index_file );
    };

    if ($@) {
        unlink($index_file);

        #$db_orth = Bio::Index::GenBank->new( -filename   => $index_file, -write_flag => 'WRITE');
        $db_orth = Bio::Index::GenBank->new( -write_flag => 'WRITE');
        $db_orth->make_index($file);
    }


    return $db_orth;
}

#
# Returns list of IDs by parsing $index->db() output. This ignores all ID
# starting with "__"!
#
sub Bio::Index::GenBank::getIDs {
    my $index = shift;

    my @keys = grep { $_ !~ /^_{2}/ } ( keys %{ $index->db() } );    # all ids

    return @keys;
}

=head2 getStat

Usage:
    
    my $arrayref = GenbankHelper::getstat($seq, $type);

Function:

  Returns assembly statistics (fraction assembled and identity) using comment
  info (written by scaffolding.pl). 


=cut

sub getStat {
    my $seq     = shift;
    my $type    = shift;

    my ($anno) = $seq->annotation()->get_Annotations('comment');
    my $string = $anno->as_text;

    my $na = "NA";
    my %hash = (
        frac_assembled => $na,
        frac_identical => $na,
        unaligned_3 => $na,
        unaligned_5 => $na,
        fragments => $na,
        bbh => $na
    );

    if ($type eq "fragments") {
        if ($string =~ /<fragments>(.+)<\/fragments>/) {
            my @e = split(", ", $1);
            $hash{"fragments"} = \@e;
        }
    } elsif ($type eq "hit-type") {
        if ($string =~ /<hit-type>(.+)<\/hit-type>/) {
            $hash{"hit-type"}  = $1;
        }
    } elsif ($type eq "best-contig") {
        if ($string =~ /<best-contig>(.+)<\/best-contig>/) {
            $hash{"best-contig"}  = $1;
        }
    } else {


    if ($string =~ /<$type>(.+)<\/$type>/) {
        my $info = $1;
        if ($info =~ /fraction assembled: (\d+\.\d+)/) {
            $hash{frac_assembled} = $1;
        }
        if ($info =~ /fraction identical: (\d+\.\d+)/) {
            $hash{frac_identical} = $1;
        }
        if ($info =~ /unaligned bases 3': (\d+)/) {
            $hash{unaligned_3} = $1;
        }
        if ($info =~ /unaligned bases 5': (\d+)/) {
            $hash{unaligned_5} = $1;
        }
    }
    }
    return \%hash;
}

1;

