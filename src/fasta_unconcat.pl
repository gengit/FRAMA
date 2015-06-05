#!/usr/bin/env perl 

=pod

=head1 INFO

Martin Bens, bensmartin@gmail.com
06/13/2014 17:00:46

=head1 DESCRIPTION

Splits single long sequence. Split based on tab-separated file with ID, start,
end in 1-based coordinates.

=head1 OPTIONS

=over 8

=item B<-file> s

Fasta file with single entry.

=item B<-index> s

Index file: ID<tab>Start<tab>End
(1-based coordinates)

=item B<-output> s

Output file. [default: STDOUT]

=back

=cut

use strict;
use warnings;

use Getopt::Long;
use File::Basename;
use Data::Dumper;
use Pod::Usage;

pod2usage(-message => "\n\tNo arguments\n", -verbose => 1) if (@ARGV == 0);

my $man  = 0;
my $help = 0;
my ($file, $index, $output);
GetOptions(
    'file|f=s'   => \$file,
    'index|i=s'  => \$index,
    'output|o=s' => \$output,
    'help|?'     => \$help,
    'man|m'      => \$man
) or pod2usage(1);

pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;

unless (-e $file && -e $index) {
    die "\n\tFasta or index file not found: Fasta:$file Index:$index\n\n";
}

my $out_fh = \*STDOUT;
if ($output) {
    open $out_fh, ">", $output or die $!;
}

my @subseqs;
open my $fh, "<", $file or die $!;
while(<$fh>) {
    chomp;
    if (/^>/) {

    } else {
        push @subseqs, $_;
    }

}
close $fh;
my $seq_string = join "", @subseqs;

open $fh, "<", $index or die $!;
while (<$fh>) {
    my ($id, $start, $end) = split "\t";
    $start--;
    my $length = $end - $start;
    toFasta($id, substr $seq_string, $start, $length);
}
close $fh;

sub toFasta {
    my ($seqName, $seq, $len) = @_;

    $len = 60 unless $len;

    print ">$seqName\n";
    while (my $chunk = substr($seq, 0, $len, "")) {
        print $chunk;
        print "\n";
    }
}

1;

__END__


