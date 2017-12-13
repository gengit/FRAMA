#!/usr/bin/env perl
#
#       AUTHOR: Martin Bens, bensmartin@gmail.com
# ORGANIZATION: FLI Jena
#

=pod

=head1 INFO

Martin Bens  (bensmartin@gmail.com)

=head1 DESCRIPTION

Find overlap between specific transcript pairs (provided by -bbh) based on
blast hits.

Output:
 1 possibleFusionID
 2 bbhID
 3 target
 4 strand
 5 possibleFusionStart
 6 possibleFusionEnd
 7 bbhStart
 8 bbhEnd
 9 Id
10 CovPosFus

=head1 OPTIONS

=over 8

=item B<-input>

Blast tabular output, as produced by blast_sort.pl. This file has to be sorted
by Target. Mapping reference transcripts to an assembly  (query: annotation,
database: assembly).

QUERY_ID     => 0
TARGET_ID    => 1
STRAND       => 2
TARGET_START => 3
TARGET_END   => 4
QUERY_COV    => 5
IDENTITY     => 6

=item B<-bbh>

Table with contigID in first column and targetID in second column. We look only
in these target transcripts for potentialy fused transcripts.

=item B<-lq>

Length of query transcripts. (first column queryID and second column
queryLength)

=item B<-max-overlap>

Maximum of allowed overlapped (in perc.) between fusion transcripts (based on
transcript size).  [int(max(bbh_length,transcript_length) * (max-overlap/100))]

=item B<-min-identity>

Minimum required identity of hit

=item B<-min-frac-size>

Minimum size of transcript which shows no overlap.

=back

=cut

use strict;
use warnings;

use Getopt::Long;
use File::Basename;
use Data::Dumper;

use Bio::Range;

use constant QUERY_ID     => 0;
use constant TARGET_ID    => 1;
use constant STRAND       => 2;
use constant TARGET_START => 3;
use constant TARGET_END   => 4;
use constant QUERY_COV    => 5;
use constant IDENTITY     => 6;

my $length_query;
my $input;
my $thr_querycov  = 100.0;
my $thr_identity  = 100.0;
my $thr_overlap   = 50.0;
my $min_frac_size = 10;
my $bbh;
my $hit_rev;
my $bestonly = 0;
my $debug = 0;

GetOptions(
    "input=s"         => \$input,
    "hit-rev=s"       => \$hit_rev,
    "bbh=s"           => \$bbh,
    "lq=s"            => \$length_query,
    "min-coverage=f"  => \$thr_querycov,
    "max-overlap=f"   => \$thr_overlap,
    "min-identity=f"  => \$thr_identity,
    "min-frac-size=i" => \$min_frac_size,
    'bestonly' => \$bestonly,
    'debug' => \$debug,
);

$thr_overlap = $thr_overlap / 100;

if (!defined $bbh || !-e $bbh) {
    die "\nNo contigs specified\n\n";
}

if (!defined $length_query || !-e $length_query) {
    die "\nProvide index file with query ID and length!\n\n";
}

# read length of targets
my %lq;
open my $fh, "<", $length_query or die "$length_query not found!";
while (<$fh>) {
    chomp;
    my ($id, $length) = split "\t";
    $lq{$id} = $length;

}
close $fh;

# hits fusion transcript on targets
my %fusion2target;
open $fh, "<", $hit_rev or die $!;
while (<$fh>) {
    chomp;
    my @e = split "\t";
    push @{$fusion2target{$e[QUERY_ID]}}, $e[TARGET_ID];
}
close $fh;

# best bidirectional hit result
my %bbh_fusion2target;
my %bbh_target2fusion;
open $fh, "<", $bbh or die $!;
while (<$fh>) {
    chomp;
    my @e = split "\t";
    # [contig] = [ortholog, startInContig, endInContig]
    $bbh_fusion2target{$e[0]} = [$e[1], $e[5], $e[6]];
    $bbh_target2fusion{$e[1]} = 1;
}
close $fh;

# HEADER
print "# possibleFusionID\tbbhID\tcontig\tstrand\tpossibleFusionStart\tpossibleFusionEnd\tbbhStart\tbbhEnd\tId\tCovPosFus\n";

my @lines;
my $bbh_line;
my $current_fusion = "";

my $inputfh;
if ($input) {
    open $inputfh, "<", $input or die $!;
} else {
    $inputfh = \*STDIN;
}
while (<$inputfh>) {
    chomp;
    my @e = split "\t";

    # new result if another target
    if ($current_fusion ne $e[TARGET_ID]) {

        # if contig has BBH, look for fusion genes
        process() if (scalar @lines > 0 && exists $bbh_fusion2target{$current_fusion});

        @lines          = ();
        $current_fusion = $e[TARGET_ID];
    }

    # collect all references that hit contig
    push @lines, [@e];
}
process() if (exists $bbh_fusion2target{$current_fusion});
close $inputfh;

sub process {

    my $bbh_id_target = $bbh_fusion2target{$current_fusion}->[0];

    my @remaining_lines = grep { $_->[QUERY_ID] ne $bbh_id_target } @lines;
    return unless (@remaining_lines);

    if ($debug) {
        print STDERR "Processing: $current_fusion-$bbh_id_target\n";
        print STDERR "\tResults: ".scalar @remaining_lines."\n";
    }

    @remaining_lines = $remaining_lines[0] if ($bestonly);

    # filter for
    # - already annotated references
    # - coverage
    # - identity
    @remaining_lines = grep {
        ($_->[QUERY_COV] >= $thr_querycov) &&
        ($_->[IDENTITY] >= $thr_identity) &&
        (!exists $bbh_target2fusion{$_->[QUERY_ID]})
    } @remaining_lines;

    print STDERR "\t\tafter filtering: ".scalar @remaining_lines."\n" if ($debug);

    return unless (@remaining_lines);

    # compute overlap always in relation to BBH range
    my $bbh_range_on_fusion = Bio::Range->new(
        -start => $bbh_fusion2target{$current_fusion}->[1],
        -end   => $bbh_fusion2target{$current_fusion}->[2]
    );

    print STDERR "\tBBH: ".$bbh_range_on_fusion->start."-".$bbh_range_on_fusion->end."\n" if ($debug);

    my $bbh_length = $lq{$bbh_id_target};
    die "Length of $bbh_id_target not found!" unless ($bbh_length);

    my @references = @{$fusion2target{$current_fusion}};
    for (@remaining_lines) {

        my $sbh_range_on_fusion = Bio::Range->new(
            -start => $_->[TARGET_START],
            -end   => $_->[TARGET_END]
        );

        print STDERR "\t\tFUSION? $_->[QUERY_ID]: ".$sbh_range_on_fusion->start."-".$sbh_range_on_fusion->end."\n" if ($debug);

        # SBH fully embedded in BBH? Skip!
        next if ($bbh_range_on_fusion->contains($sbh_range_on_fusion));

        # calculate unique regions in BBH, region in common und unique region of SBH
        my ($bbh_unique, $common, $range_unique) = $bbh_range_on_fusion->overlap_extend($sbh_range_on_fusion);

        # skip lines below minimum fragment size
        next if ($range_unique < $min_frac_size);

        # we might have a better ortholog here than the annotated BBH. Skip.
        next if ($bbh_unique == 0);

        my $sbh_length = $lq{$_->[QUERY_ID]};
        die "Length of $_->[QUERY_ID] not found!" unless ($sbh_length);

        my $allowed_overlap = int(($thr_overlap * getMax($bbh_length, $sbh_length)) + 0.5);

        if ($common <= $allowed_overlap) {
            my $found = 0;
            if (@references) {
                for my $ref (@references) {
                    if ($ref eq $_->[QUERY_ID]) {
                        $found = 1;
                        last;
                    }
                }
            }

            if ($found) {
                print join("\t", $_->[QUERY_ID], $bbh_id_target, $_->[TARGET_ID], $_->[STRAND], $_->[TARGET_START], $_->[TARGET_END], $bbh_range_on_fusion->start, $bbh_range_on_fusion->end, $_->[IDENTITY], $_->[QUERY_COV]);
                print "\n";
                # ignore hits to additional contigs
                $bbh_target2fusion{$_->[QUERY_ID]} = 1;
            }
        }
    }
}

sub Bio::Range::overlap_extend {
    my ($a, $b) = @_;

    $a->throw("start is undefined")       unless defined $a->start;
    $a->throw("end is undefined")         unless defined $a->end;
    $b->throw("Not a Bio::RangeI object") unless $b->isa('Bio::RangeI');
    $b->throw("start is undefined")       unless defined $b->start;
    $b->throw("end is undefined")         unless defined $b->end;

    if (!$a->overlaps($b)) {
        return ($a->length, 0, $b->length);
    }

    my ($au, $bu) = (0, 0);
    if ($a->start < $b->start) {
        $au = $b->start - $a->start;
    } else {
        $bu = $a->start - $b->start;
    }

    if ($a->end > $b->end) {
        $au += $a->end - $b->end;
    } else {
        $bu += $b->end - $a->end;
    }

    my $intersect = $a->intersection($b);
    if (!$intersect) {
        warn("no intersection\n");
        return ($au, 0, $bu);
    } else {
        my $ie = $intersect->end;
        my $is = $intersect->start;
        return ($au, $ie - $is + 1, $bu);
    }
}

sub getMax {
    my ($x, $y) = @_;

    if ($x > $y) {
        return $x
    } else {
        return $y
    }
}

1;

__END__

