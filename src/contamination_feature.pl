#!/usr/bin/env perl

=pod

=head1 INFO

Martin Bens, bensmartin@gmail.com
01/04/14 16:11:48

=head1 DESCRIPTION

Adds misc_feature for regions which might originate from vector sequences.

=head1 OPTIONS

=over 8

=item B<-input>

Input file.

=item B<-output>

Target file.

=item B<-contamination>

CSV as produced by vecscreen.pl

=back

=cut

use strict;
use warnings;

use Getopt::Long;
use File::Basename;
use Data::Dumper;
use Pod::Usage;

use Bio::SeqIO;

pod2usage(-message => "\n\tNo arguments\n", -verbose => 1) if (@ARGV == 0);

my $man  = 0;
my $help = 0;
my ($input, $output, $contamination, $reindex);
$reindex = 0;
GetOptions(
    'input|i=s'         => \$input,
    'output|o=s'        => \$output,
    'contamination|c=s' => \$contamination,
    'reindex'           => \$reindex,
    'help|?'            => \$help,
    'man|m'             => \$man
) or pod2usage(1);

pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;

my %contamination;
my $current_id;
open my $fh, "<", $contamination or die $!;
while (<$fh>) {
    chomp;
    my @e = split "\t";
    push @{$contamination{$e[0]}}, [@e];
}
close $fh;

if ((keys %contamination) == 0) {
    system("cp $input $output");
} else {
    my $genbank = Bio::SeqIO->new(-file => $input,        -format => 'genbank');
    my $out     = Bio::SeqIO->new(-file => ">" . $output, -format => 'genbank');
    while (my $seq = $genbank->next_seq) {
        if (exists $contamination{$seq->id}) {

            my @contaminations = @{$contamination{$seq->id}};

            for my $e (@contaminations) {
                my $id = $e->[4];
                if ($e->[4] =~ /.*\|(.+?)\:/) {
                    $id = $1;
                }
                my $feature = Bio::SeqFeature::Generic->new(
                    -start   => $e->[6],
                    -end     => $e->[7],
                    -primary => 'misc_feature',
                    -tag     => {
                        match       => $id,
                        description => $e->[5],
                        note        => "$e->[2]:$e->[3]:$e->[8]"
                    }
                );
                $seq->add_SeqFeature($feature);
            }
        }

        $out->write_seq($seq);
    }
}

1;

__END__

