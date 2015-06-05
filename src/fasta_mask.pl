#!/usr/bin/env perl

use warnings;
use strict;

=pod

=head1 INFO

Martin Bens, bensmartin@gmail.com
2014-08-06

=head1 DESCRIPTION

Little helper if RepeatMasker finished and you still have the intermediate file
containing masking region, but not the actual masked fasta file (or it was not
written). Extract column with name, start and end and specify the resulting
table with "masking" parameter.

=head1 OPTIONS

=over 8

=item B<-masking>

Table with SequenceID, Start and end.

=back

=cut

use Pod::Usage;
use Getopt::Long;
pod2usage(-message => "\n\tNo arguments\n", -verbose => 1) if (@ARGV == 0);

my $man  = 0;
my $help = 0;
my $region_file;
GetOptions(
    'masking=s' => \$region_file,
    'help|?' => \$help,
    'man|m'  => \$man
) or pod2usage(1);

pod2usage(1) if $help;
pod2usage( -verbose => 2 ) if $man;

my %masking;
open my $fh, "<", $region_file or die $!;
while(<$fh>) {
    chomp;
    my @e = split;
    next unless (@e == 3 && $e[1] =~ /^\d/ && $e[2] =~ /^\d/);
    push @{$masking{$e[0]}}, [$e[1]-1, $e[2] - $e[1] + 1];
}
close $fh;

my $current_name;
my @current_sequence;

while (<STDIN>) {
    chomp;
    if (/^>(.+?)(\s|$)/) {
        my $string = join "", @current_sequence if (@current_sequence);
        @current_sequence = ();
        if (defined $current_name && $masking{$current_name}) {
            for my $region (@{$masking{$current_name}}) {
                substr $string, $region->[0], $region->[1], lc(substr $string, $region->[0], $region->[1]);
            }
        }
        toFasta($current_name, $string, 60) if ($current_name);
        $current_name = $1;

    } else {
        push @current_sequence, $_;
    }
}
my $string = join "", @current_sequence;
if (defined $current_name && $masking{$current_name}) {
    for my $region (@{$masking{$current_name}}) {
        substr $string, $region->[0], $region->[1], lc(substr $string, $region->[0], $region->[1]);
    }
}
toFasta($current_name, $string, 60);

sub toFasta {
    my ($seqName, $seq, $len) = @_;

    $len = 60 unless $len;

    print ">$seqName\n";
    while (my $chunk = substr($seq, 0, $len, "")) {
        print $chunk."\n";
    }
}

