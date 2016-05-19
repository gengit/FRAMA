#!/usr/bin/env perl

=pod

=head1 INFO

Martin Bens, bensmartin@gmail.com
03/19/2014 09:45:53

=head1 DESCRIPTION

Combines all HSP of one hit. Total score is sum of all HSP-scores.
Additionally, highest HSP-score and covered regions on query and target are
added. Identity of hit is mean of all HSPs.

Format: 
 1 QueryID,
 2 TargetID,
 3 #HSPs,
 4 #Evalue,
 5 #TotalScore,
 6 #HighScore,
 7 #Group,
 8 #Strand,
 9 QueryMinStart,
 10 QueryMaxEnd,
 11 TargetMinStart,
 12 TargetMaxEnd,
 13 CovBpQuerp,
 14 CovBpTarget,
 15 CovPercQuery,
 16 CovPercTarget,
 17 MeanIdentity

=head1 OPTIONS

=over 8

=item B<-input>

WU-Blast result in table format with topcombo and links column.

=item B<-lq> 

File with length of queries (Column1: ID, Column2: length)

=item B<-lt> 

File with length of target (Column1: ID, Column2: length)

=item B<-output> 

Output file (default: STDOUT)

=back

=cut

use strict;
use warnings;

use FindBin;
use lib "$FindBin::Bin/../lib/perl";

use Getopt::Long;
use Data::Dumper;
use Pod::Usage;
use Set::IntSpan;

use constant QUERY_ID     => 0;
use constant TARGET_ID    => 1;
use constant EVALUE       => 2;
use constant SCORE        => 4;
use constant STRAND       => 16;
use constant QUERY_START  => 17;
use constant QUERY_END    => 18;
use constant TARGET_START => 20;
use constant TARGET_END   => 21;
use constant GROUP        => 22;
use constant IDENTITY     => 10;

pod2usage(-message => "\n\tNo arguments.\n", -verbose => 1) if (@ARGV == 0);

my $man  = 0;
my $help = 0;

my ($file_blast, $file_query, $file_target, $output);
GetOptions(
    'input|i=s'          => \$file_blast,
    'length-query|lq=s'  => \$file_query,
    'length-target|lt=s' => \$file_target,
    'output|o=s'         => \$output,
    'help|?'             => \$help,
    'man|m'              => \$man
) or pod2usage(1);

pod2usage(1) if $help;
pod2usage(-verbose => 2 ) if $man;

unless (-e $file_query && -e $file_target && -e $file_blast) {
    die "At least one of the following files was not found : \n\t$file_query,\n\t$file_target,\n\t$file_blast";
}
if (-z $file_query || -z $file_target || -z $file_blast) {
    die "At least one of the following files is empty: \n\t$file_query,\n\t$file_target,\n\t$file_blast";
}

my (%lq, %lt);
open my $fh, "<", $file_query or die "Could not open file: $!";
while (<$fh>) {
    chomp;
    my ($id, $length) = split "\t";
    $lq{$id} = $length;
}
close $fh;
open $fh, "<", $file_target or die "Could not open file: $!";
while (<$fh>) {
    chomp;
    my ($id, $length) = split "\t";
    $lt{$id} = $length;
}
close $fh;

my $outfh;
if ($output) {
    open $outfh, ">", $output or die $!;
} else {
    $outfh = \*STDOUT;
}

my ($c_query, $c_target, $c_group) = ("", "", "");
my @lines;
open $fh, "<", $file_blast or die "Could not open file: $!";
while (<$fh>) {
    chomp;
    my @e = split "\t";
    unless ($c_group eq $e[GROUP] && $c_target eq $e[TARGET_ID] && $c_query eq $e[QUERY_ID]) {
        process() if (@lines > 0);
        @lines    = ();
        ($c_query, $c_target, $c_group) = @e[QUERY_ID, TARGET_ID, GROUP];
    }
    push @lines, [@e[EVALUE,SCORE,STRAND,QUERY_START,QUERY_END,TARGET_START,TARGET_END,IDENTITY]];
}
process();
close $fh;
close $outfh;

sub process {

    my $hsps      = scalar @lines;
    my $firstline = shift @lines;

    my $total_score = $firstline->[1];
    my $high_score  = $firstline->[1];
    my $identity    = $firstline->[7];

    my $query_l  = $lq{$c_query};
    my $target_l = $lt{$c_target};

    unless ($query_l) {
        print STDERR "Length of query not found: $c_query\n";
        return;
    }
    unless ($target_l) {
        print STDERR "Length of target not found: $c_target\n";
        return;
    }

    my ($index_start, $index_end) = (3, 4);

    # colinear hsps => all same strand
    if ($firstline->[2] eq "-1") {
        ($index_start, $index_end) = ($index_end, $index_start);
    }

    my $query_mapped_l;
    my $target_mapped_l;

    my $query_start;
    my $query_end;

    my $target_start;
    my $target_end;

    if (@lines > 0) {
        my @query_ranges  = ($firstline->[$index_start] .. $firstline->[$index_end]);
        my @target_ranges = ($firstline->[5] .. $firstline->[6]);

        for (@lines) {
            $total_score += $_->[1];
            $high_score = ($high_score < $_->[1]) ? $_->[1] : $high_score;
            push @query_ranges,  ($_->[$index_start] .. $_->[$index_end]);
            push @target_ranges, ($_->[5] .. $_->[6]);
            $identity += $_->[7];
        }
        my $span_query  = Set::IntSpan->new(@query_ranges);
        my $span_target = Set::IntSpan->new(@target_ranges);

        $query_mapped_l  = $span_query->size;
        $target_mapped_l = $span_target->size;

        $query_start = $span_query->min;
        $query_end   = $span_query->max;

        $target_start = $span_target->min;
        $target_end   = $span_target->max;
    } else {
        $query_start = $firstline->[$index_start];
        $query_end   = $firstline->[$index_end];

        $target_start = $firstline->[5];
        $target_end   = $firstline->[6];

        $query_mapped_l  = $query_end - $query_start + 1;
        $target_mapped_l = $target_end - $target_start + 1;
    }

    print $outfh join(
        "\t",
        $c_query,
        $c_target,
        $hsps,              # hsps
        $firstline->[0],    # evalue
        $total_score,
        $high_score,
        $c_group,
        $firstline->[2],    #strand
        $query_start,
        $query_end,
        $target_start,
        $target_end,
        $query_mapped_l,
        $target_mapped_l,
        sprintf("%.2f", (($query_mapped_l / $query_l) * 100)),      # perc. covered query
        sprintf("%.2f", (($target_mapped_l / $target_l) * 100)),    # perc. covered target
        sprintf("%.2f", ($identity / $hsps))                        # mean identity
    );
    print $outfh "\n";
}

1;

__END__

