#!/usr/bin/env perl

=pod

=head1 INFO

Martin Bens, bensmartin@gmail.com
12/18/2013 17:24:53

=head1 DESCRIPTION

=head1 OPTIONS

=over 8

=item B<-input>

Input fasta file.

=item B<-genbank>

Related genbank file

=item B<-output>

Output file (fasta).

=back

=cut

use strict;
use warnings;

use Getopt::Long;
use File::Basename;
use Data::Dumper;
use Pod::Usage;

use Cwd qw(realpath);
BEGIN {
    my ($mypath) = realpath($0)=~m/(^.*)\//;
    push @INC, "$mypath/../lib/perl";
}

use GenbankHelper;

pod2usage(-message => "\n\tNo arguments\n", -verbose => 1) if (@ARGV == 0);

my $man  = 0;
my $help = 0;
my ($input, $genbank_file, $output);
GetOptions(
    'input=s'   => \$input,
    'genbank=s' => \$genbank_file,
    'output=s'  => \$output,
    'help|?'    => \$help,
    'man|m'     => \$man
) or pod2usage(1);

pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;

if (!-e $input && -z $input) {
    die "Input not found: $input\n";
}
if (!-e $genbank_file && -z $genbank_file) {
    die "Input not found: $genbank_file\n";
}

my $genbank = GenbankHelper::getIndex($genbank_file);

my $seq_io  = Bio::SeqIO->new(-file => $input,        -format => 'fasta');
my $seq_out = Bio::SeqIO->new(-file => ">" . $output, -format => 'fasta');

while (my $seq = $seq_io->next_seq) {
    my $gb_seq = $genbank->get_Seq_by_id($seq->id);

    my ($cds) = $gb_seq->get_SeqFeatures('CDS');

    if ($cds) {
        $seq_out->write_seq($seq->trunc($cds->start, $cds->end));
    }
}

1;

__END__

