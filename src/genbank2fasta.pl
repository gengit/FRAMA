#!/usr/bin/env perl
#
#       AUTHOR: Martin Bens, bensmartin@gmail.com
# ORGANIZATION: FLI Jena
#      CREATED: 02/21/13 15:15:06
#

use strict;
use warnings;

use FindBin;
use lib "$FindBin::Bin/../lib/perl";

use Getopt::Long;
use Data::Dumper;
use Bio::SeqIO;
use Pod::Usage;

use GenbankHelper;

my $input;
my $output;
my $ids;
my $id;
my $format = "fasta";
my $i_name = "acc";
my $o_name = "acc";
my $feature;
my $tag;
my $value;
my $help             = 0;
my $man              = 0;
my $feature_sequence = 0;
my $complete         = 0;
my $protein          = 0;
my $molecule         = 0;

if (@ARGV == 0) {
    pod2usage();
    exit;
}

GetOptions(
    "input|i=s"           => \$input,
    "output|o=s"          => \$output,
    "feature|feat=s"      => \$feature,
    "feature-sequence|fs" => \$feature_sequence,
    "molecule=s"          => \$molecule,
    "format|f=s"          => \$format,
    "ids=s"               => \$ids,
    "id=s"                => \$id,
    "o-name|on=s"         => \$o_name,
    "i-name|in=s"         => \$i_name,
    "tag|t=s"             => \$tag,
    "value|v=s"           => \$value,
    "complete|c"          => \$complete,
    "protein|p"           => \$protein,
    "help|?"              => \$help,
    "man|?"               => \$man,
) or die;

pod2usage() if ($help);
pod2usage(-verbose => 2) if ($man);

die "\n\tInput file not found: $input\n\n" unless (-e $input);

my $io = Bio::SeqIO->new(-file => $input, -format => 'GenBank');
my $builder = $io->sequence_builder();
$builder->want_none();
$builder->add_wanted_slot('seq', 'features');

my $out;
if ($output) {
    $out = Bio::SeqIO->new(-file => ">$output", -format => $format, -flush => 0);
} else {
    $out = Bio::SeqIO->new(-fh => \*STDOUT, -format => $format, -flush => 0);
}

if ($protein) {
    $feature          = 0;
    $feature_sequence = 1;
}

# we are going to change transcript accession in order to allow different
# identifier in fasta header
$out->preferred_id_type('accession');

my %ids;
if (defined $id) {
    $ids{$id} = 1;
} elsif (defined $ids) {
    open my $fh, "<", $ids or die $!;
    while (<$fh>) {
        chomp;
        $ids{$_} = 1;
    }
    close $fh;
}

my $output_name = 0;
if ($o_name eq "acc") {
    $output_name = 0;
} elsif ($o_name eq "acc.ver") {
    $output_name = 1;
} elsif ($o_name eq "symbol") {
    $output_name = 2;
} elsif ($o_name eq "pipeline_acc.ver") {

    # to extract RefSeqID from pipeline genbank
    $output_name = 3;
} elsif ($o_name eq "pipeline_acc") {

    # to extract RefSeqID from pipeline genbank
    $output_name = 4;
}

my $input_name = 0;
if ($i_name eq "acc") {
    $input_name = 0;
} elsif ($i_name eq "acc.ver") {
    $input_name = 1;
} elsif ($i_name eq "symbol") {
    $input_name = 2;
}

while (my $seqobj = $io->next_seq) {
    my $seqid = $seqobj->accession;
    $seqid = $seqobj->id if ($seqid eq "unknown"); 
    my $outseq;

    if ($molecule) {
        next if ($seqobj->molecule ne $molecule);
    }
    if ($feature) {
        my @feat = $seqobj->get_SeqFeatures($feature);
        next unless (@feat > 0);

        my $location = $feat[0]->location;
        if ($complete) {
            next
              unless ($location->start_pos_type eq "EXACT"
                && $location->start_pos_type eq "EXACT");
        }

        if (defined $tag) {
            for my $feat (@feat) {
                if ($feat->has_tag($tag)) {
                    if (defined $value) {
                        my @tagval = $feat->get_tag_values($tag);
                        for my $tagvalue (@tagval) {
                            if ($tagvalue =~ /$value/) {
                                if ($feature_sequence) {
                                    $outseq = $feat->seq;
                                } else {
                                    $outseq = $seqobj;
                                }
                            }
                        }
                    } else {
                        if ($feature_sequence) {
                            $outseq = $feat->seq;
                        } else {
                            $outseq = $seqobj;
                        }
                    }
                }
            }
        } else {
            if ($feature_sequence) {
                $outseq = $feat[0]->spliced_seq;
            } else {
                $outseq = $seqobj;
            }
        }
    } elsif ($protein) {
        my @feat = $seqobj->get_SeqFeatures("CDS");
        next unless (@feat > 0);
        if ($feat[0]->has_tag("translation")) {
            my ($string) = $feat[0]->get_tag_values('translation');
            $outseq = Bio::Seq->new(-id => $seqobj->id, -seq => $string);
        }
    } else {
        $outseq = $seqobj;
    }

    next unless (defined $outseq);

    if ($input_name) {
        if ($input_name == 1) {
            $seqid = $seqid . "." . $seqobj->version;
        } elsif ($input_name == 2) {
            $seqid = GenbankHelper::getSymbol($seqobj);
            $seqid = $seqobj->accession unless ($seqid);
            $seqid = $seqobj->id if ($seqid eq "unknown");
        }
    }

    next if (%ids && !exists $ids{$seqid});

    if ($output_name == $input_name) {
        $outseq->accession_number($seqid);
    } else {
        if ($output_name == 1) {
            $seqid = $seqid . "." . $seqobj->version;
        } elsif ($output_name == 2) {
            $seqid = GenbankHelper::getSymbol($seqobj);
            $seqid = $seqobj->accession unless ($seqid);
            $seqid = $seqobj->id if ($seqid eq "unknown");
        } elsif ($output_name == 3) {
            my ($feature) = $seqobj->get_SeqFeatures('gene');
            if ($feature->has_tag('inference')) {
                my ($value) = $feature->get_tag_values('inference');
                my @e = split ":", $value;
                $seqid = $e[4];
            }
        } elsif ($output_name == 4) {
            my ($feature) = $seqobj->get_SeqFeatures('gene');
            if ($feature->has_tag('inference')) {
                my ($value) = $feature->get_tag_values('inference');
                my @e = split ":", $value;
                ($seqid) = split(/\./, $e[4]);
            }

        }
        $outseq->accession_number($seqid);
    }

    $out->write_seq($outseq);
}

1;

__END__

=head1 DESCRIPTION
    
Converts genbank file to fasta format.

=head1 SYNOPSIS

Convert genbank to fasta.
        
    genbank2fasta.pl -i input.gbk -o output.fa

Convert genbank to fasta, but use symbol as header.
        
    genbank2fasta.pl -i input.gbk -o output.fa -o-name symbol

Extract CDS sequence and write to fasta file.

    genbank2fasta.pl -i input.gbk -o output.fa -format fasta -feature CDS -feature-sequence

Extract only sequences with CDS annotated and write to fasta file.

    genbank2fasta.pl -i input.gbk -o output.fa -format fasta -feature CDS 

Extract only sequence if miRBase-ID.

    genbank2fasta.pl -i input.gbk -o output.fa -format fasta -feature gene -tag db_xref -value miRBase


=head1 OPTIONS

=over 8

=item B<-input|-i> input.gbk
    
Genbank file.

=item B<-output|-o> output.fa

Output file.

=item B<-o-name|-on> acc|acc.ver|symbol

default: acc

    acc      ACCESSION as header
             (e.g. ">NM_000000 description")
    acc.ver  ACCESSION.VERSION as header
             (e.g. ">NM_00000.1 description")
    symbol   gene symbol as annotated in "gene" feature.
                 (e.g. ">A2M description")

=item B<-i-name|-in> acc|acc.ver|symbol

default: acc

ID type of input IDs (see -o-name).

=item B<-feature|-feat> CDS|gene|...

Entry must include specified feature in order to be converted to fasta.

=item B<-feature-sequence|-fs> 

Truncate sequence by feature coordinates (see -feature). Useful to extract CDS
region for instance.

=item B<-ids> file.txt
    
IDs to extract from input file. One ID per line (see -i-name).

=item B<-id> A2M

ID to extract from input (see -i-name).

=item B<-tag|-t> db_xref|pseudo|...

Only convert entry if feature has specified tag.

=item B<-value|-v> GeneID|...

Only convert entry if value of -tag matches specified text (regexp).

=item B<-format> fasta|genbank|...

Bioperl compatible sequence format. [default: fasta]

