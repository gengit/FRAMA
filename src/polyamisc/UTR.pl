#!/usr/bin/env perl 

=pod

=head1 INFO

Martin Bens, bensmartin@gmail.com
07/10/2014 17:19:24

=head1 DESCRIPTION

Computes identity of UTR aligned to contig (first sequence must be UTR).
Identity is computed only for aligned region of UTR. Contribution of end-gaps
in contig to identity can be specified by -end-gap.
=head1 OPTIONS

=over 8

=item B<-input>

UTR alignment

=item B<-output>

Output file. [default: STDOUT]

=item B<-end-gap> 

Reward for end-gaps. 

0 => count as mismatch
1 => count as match

[default: 0.5]

=back

=cut

use strict;
use warnings;

use Getopt::Long;
use File::Basename;
use Data::Dumper;
use Pod::Usage;
use Bio::AlignIO;

pod2usage(-message => "\n\tNo arguments\n", -verbose => 1) if (@ARGV == 0);

my $man  = 0;
my $help = 0;
my ($input, $output);
my $utr     = "UTR";
my $format  = "fasta";
my $end_gap = 0.5;

GetOptions(
    'input|i=s'    => \$input,
    'output|o=s'   => \$output,
    'end-gap|eg=f' => \$end_gap,
    'format=s'     => \$format,
    'help|?'       => \$help,
    'man|m'        => \$man
) or pod2usage(1);

pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;

die "Alignment file not found: $input" unless ($input && -e $input);
my $aln = Bio::AlignIO->new(-file => $input, -format => $format)->next_aln;
die "Alignment not readable" unless ($aln);

my $utr_seq    = $aln->get_seq_by_pos(1);
my $contig_seq = $aln->get_seq_by_pos(2);

die "Not enough sequences" unless ($utr_seq && $contig_seq);

my $utr_id    = $utr_seq->id;
my $contig_id = $contig_seq->id;

my ($left_overhang, $right_overhang, $start, $end);
$left_overhang  = 0;
$right_overhang = 0;
$start          = 1;
$end            = $aln->length;

if ($utr_seq->seq =~ /^(-+)/) {
    $left_overhang = length($1);
    $start         = length($1) + 1;
}
if ($utr_seq->seq =~ /(-+)$/) {
    $right_overhang = length($1);
    $end -= length($1);
}

# remove end gaps
if ($left_overhang != 0 || $right_overhang != 0) {
    $aln->verbose(-1);    # ignore warning
    $aln = $aln->slice($start, $end);
}

my $end_pos  = 0;
my $identity = 0;

# possible that UTR does not align with contig
if ($contig_seq->location_from_column($end)) {
    $end_pos  = $contig_seq->location_from_column($end)->start;
    $identity = identity($aln);
} 

my $fh;
if ($output) {
    open $fh, ">", $output or die $!;
} else {
    $fh = \*STDOUT;
}

print $fh "###\n";
print $fh "# sequence $contig_id\n";
print $fh "#\n";
printf $fh "%s\t%.3f\n", $end_pos, $identity;

close $fh;

sub identity {
    my ($utr, $contig) = $aln->each_seq;
    my $length = $aln->length;

    my $end          = $length;
    my $num_end_gaps = 0;
    if ($contig->seq =~ /(-+)$/) {
        $end          = $length - length($1);
        $num_end_gaps = length($1);
    }

    my @utrChar    = split "", uc($utr->seq);
    my @contigChar = split "", uc($contig->seq);
    my $m;
    for (my $i = 0; $i < $end; $i++) {
        if ($utrChar[$i] eq $contigChar[$i]) {
            $m++;
        }
    }
    return (($m + ($end_gap * $num_end_gaps)) / $length) * 100;
}

1;

__END__

