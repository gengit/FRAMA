#!/usr/bin/env perl

=pod

=head1 INFO

Martin Bens, bensmartin@gmail.com
01/31/13 16:52:23

=head1 DESCRIPTION

Extracts line with first occurence of a value in a specified column. This can
therefore be used to extract best hits, for instance for BLAT and Blast, which
should previously be sorted by your favoured best hit criteria.

=head1 OPTIONS

=over 8

=item B<-input> s

Input file. (tab separated) [default: STDIN]

=item B<-output> s

Output file. [default: STDOUT]

=item B<-column> i

Specify column. [default: 1]

=item B<-multiple> i

Keep first "i" occurences.


=back

=cut

use strict;
use warnings;

use Getopt::Long;
use File::Basename;
use Data::Dumper;
use Pod::Usage;

my ($input, $output, $column);

$column = 1;

my $man   = 0;
my $help  = 0;
my $multi = 0;
GetOptions(
    "input=s"    => \$input,
    "output=s"   => \$output,
    "column=i"   => \$column,
    'help|?'     => \$help,
    'multiple=i' => \$multi,
    man          => \$man
);

pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;

# 1-based to 0-based
$column--;

my %queries;
my $inputfh;
my $outfh;

if ($input) {
    die "ERROR: File not found: $input\n" unless (-e $input);
    open $inputfh, "<", $input or die $!;
} else {
    $inputfh = \*STDIN;
}
if ($output) {
    open $outfh, ">", $output or die $!;
} else {
    $outfh = \*STDOUT;
}
while (<$inputfh>) {
    my @e        = split;
    my $query_id = $e[$column];

    if ($multi) {
        next if (exists $queries{$query_id} && $queries{$query_id} == $multi);
    } else {
        next if (exists $queries{$query_id});
    }
    $queries{$query_id}++;

    print $outfh $_;
}
close $outfh;
close $inputfh;

1;

__END__


