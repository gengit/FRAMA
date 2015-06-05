#!/usr/bin/env perl 

=pod

=head1 INFO

Martin Bens, bensmartin@gmail.com
06/26/2014 16:35:36

=head1 DESCRIPTION

1-based !

This script checks soft clipped regions for poly(A)-tails (mininum number of As
can be specified). You can further specify the mininum ratio of clipped As in
contrast all clipped bases.

=head1 OPTIONS

=over 8

=item B<-i>

BAM/SAM file with locally aligned reads. 

=item B<-ratio> f

default: 0.8 (= at least 80% As)

=item B<-min-A> i

Minimum number of As in clipped region.

=item B<-seq_id> 

(optional) print results only for specific seq_id

=item B<-revcom> 

(optional) look for poly(A)-heads instead.

=item B<-verbose>

(optional) prints details for all clipped reads

=back

=cut

use strict;
use warnings;

use Getopt::Long;
use File::Basename;
use Data::Dumper;
use Pod::Usage;

pod2usage(-message => "\n\tNo arguments\n", -verbose => 1) if (@ARGV == 0);

use constant CIGAR  => 5;
use constant SEQ    => 9;
use constant FLAG   => 1;
use constant TARGET => 2;
use constant POS    => 3;

my $man  = 0;
my $help = 0;
my ($file_sam, $seq_id);
my $ratio   = 0.8;
my $min_A   = 2;
my $verbose = 0;
my $revcom = 0;
GetOptions(
    'input|i=s' => \$file_sam,
    'seq_id=s'  => \$seq_id,
    'min-A=i'   => \$min_A,
    'ratio=f'   => \$ratio,
    'verbose'   => \$verbose,
    'revcom' => \$revcom,
    'help|?'    => \$help,
    'man|m'     => \$man
) or pod2usage(1);

pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;

unless (defined($file_sam) and -e $file_sam) {
    die "BAM file not found: $file_sam\n";
}

my %result;
my $fh;
if ($file_sam =~ /\.sam$/) {
    open $fh, "<", $file_sam or die $!;
} else {
    open $fh, "-|", "samtools view $file_sam" or die $!;
}
if ($revcom) {
    while(<$fh>) {
        next if (/^@/);

        my @e = split;
        next if ($seq_id && $e[TARGET] ne $seq_id);
        $result{$e[TARGET]} ||= {};

        # read sequence is already orientated
        # my $flag_bin = sprintf "%b", $e[FLAG];

        if ($e[CIGAR] =~ /(^\d{$min_A,})S/) {
            my $clipped       = $1;
            my $clipped_bases = substr $e[SEQ], $clipped;
            my $num_A         = $clipped_bases =~ tr/T//;
            my $end           = 0;
            if ($num_A > $min_A && ($num_A / $clipped) >= $ratio) {
                print STDERR "T\t" if ($verbose);
                my $length = 0;
                while ($e[CIGAR] =~ /(\d+)(\w)/g) {

                    # matches and deletion account for end position
                    if ($2 eq "M" || $2 eq "D") {
                        $length += $1;
                    }
                }

                # last aligned base
                $end = ($e[POS] - $clipped);
                $result{$e[TARGET]}->{$end}++;
            } else {
                print STDERR "F\t" if ($verbose);
            }
            print STDERR $e[TARGET] . "\t$num_A\t" . ($end) . "\t$clipped_bases\t$clipped\n" if ($verbose);
        }

    }
} else {


    while (<$fh>) {
        next if (/^@/);

        my @e = split;
        next if ($seq_id && $e[TARGET] ne $seq_id);
        $result{$e[TARGET]} ||= {};

        # read sequence is already orientated
        # my $flag_bin = sprintf "%b", $e[FLAG];

        if ($e[CIGAR] =~ /(\d{$min_A,})S$/) {
            my $clipped       = $1;
            my $clipped_bases = substr $e[SEQ], -$clipped;
            my $num_A         = $clipped_bases =~ tr/A//;
            my $end           = 0;
            if ($num_A > $min_A && ($num_A / $clipped) >= $ratio) {
                print STDERR "T\t" if ($verbose);
                my $length = 0;
                while ($e[CIGAR] =~ /(\d+)(\w)/g) {

                    # matches and deletion account for end position
                    if ($2 eq "M" || $2 eq "D") {
                        $length += $1;
                    }
                }

                # last aligned base
                $end = ($e[POS] + $length - 1);
                $result{$e[TARGET]}->{$end}++;
            } else {
                print STDERR "F\t" if ($verbose);
            }
            print STDERR $e[TARGET] . "\t$num_A\t" . ($end) . "\t$clipped_bases\t$clipped\n" if ($verbose);
        }

    }
    close $fh;
}

for my $seq_id (keys %result) {
    my %positions = %{$result{$seq_id}};

    print "###\n";
    print "# sequence $seq_id\n";
    print "#\n";

    for (sort { $a <=> $b } keys %positions) {
        print join("\t", $_, $positions{$_}) . "\n";
    }
}

1;

__END__


USING BIO::DB::SAM

    my $sam = Bio::DB::Sam->new(-bam => $file_sam);

    my @alignments = $sam->get_features_by_location(-seq_id => $seq_id);

    print STDERR "polyA\tseqID\tnumA\tstartClip\tclipped_bases\tlength\n";
    for my $a (@alignments) {

        if ($a->cigar_str =~ /(\d{$min_A,})S$/) {
            my $clipped       = $1;
            my $clipped_bases = substr $a->query->dna, -$clipped;
            my $num_A         = $clipped_bases =~ tr/A//;
            if ($num_A > $min_A && ($num_A / $clipped) >= $ratio) {
                print STDERR "T\t" if ($verbose);
                $result{$a->seq_id}->{$a->end}++;
            } else {
                print STDERR "F\t" if ($verbose);
            }
            print STDERR $a->seq_id."\t$num_A\t".($a->end)."\t$clipped_bases\t$clipped\n" if ($verbose);
        }
    }
}
