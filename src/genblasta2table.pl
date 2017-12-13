#!/usr/bin/env perl

=pod

=head1 INFO

Martin Bens, bensmartin@gmail.com
07/15/2014 18:06:35

=head1 DESCRIPTION

=head1 OPTIONS

=over 8

=item B<-input>

genblast result

=item B<-output>

outputs first rank only

=back

=cut

use strict;
use warnings;

use Cwd qw(realpath);
BEGIN {
    my ($mypath) = realpath($0)=~m/(^.*)\//;
    push @INC, "$mypath/../lib/perl";
}

use Getopt::Long;
use File::Basename;
use Data::Dumper;
use Pod::Usage;

pod2usage( -message => "\n\tNo arguments\n", -verbose => 1 ) if ( @ARGV == 0 );

my $man  = 0;
my $help = 0;
my ( $input, $output );
my $header      = 0;
my $score_strat = "half";
my ( $lq, $lt );
GetOptions(
    'input|i=s'  => \$input,
    'output|o=s' => \$output,
    'score=s'    => \$score_strat,
    'header|h!'  => \$header,
    'lq=s'       => \$lq,
    'lt=s'       => \$lt,
    'help|?'     => \$help,
    'man|m'      => \$man
) or pod2usage(1);

pod2usage(1) if $help;
pod2usage( -verbose => 2 ) if $man;

sub getMin {
    my $min = shift;

    for (@_) {
        $min = $min < $_ ? $min : $_;
    }
    return $min;
}

sub getMax {
    my $max = shift;

    for (@_) {
        $max = $max > $_ ? $max : $_;
    }
    return $max;
}

unless ( -e $lq && -e $lt ) {
    die "File with length of query or target missing\n";
}

my %length_query;
open my $fh, "<", $lq or die "Can't open file for reading: $lq\n";
while (<$fh>) {
    chomp;
    my ( $query, $length ) = split "\t";
    $length_query{$query} = $length;
}
close $fh;

my %length_target;
open $fh, "<", $lt or die "Can't open file for reading: $lq\n";
while (<$fh>) {
    chomp;
    my ( $query, $length ) = split "\t";
    $length_target{$query} = $length;
}
close $fh;

my $out;
if ($output) {
    open $out, ">", $output or die $!;
} else {
    $out = \*STDOUT;
}

if ($header) {
    print $out join "\t",
      qw/query target hsps evalue totalscore highscore group strand queryminstart querymaxend targetminstart targetmaxend covBpQ covBpT covPercQ covPercT meanIdent/;
    print $out "\n";
}

open $fh, "<", $input or die $!;
while (<$fh>) {
    chomp;
    my @e = split "\\|";

    if ( @e == 6 ) {
        my ( $target, $target_se ) = split ":", $e[1];

        my ( $target_s, $target_e ) = ( 0, 0 );
        if ( $target_se =~ /(\d+)\.\.(\d+)/ ) {
            $target_s = $1;
            $target_e = $2;
        }

        my $query_cover      = -1;
        my $query_cover_perc = -1;
        if ( $e[3] =~ /gene cover:(\d+)\((.+)\%\)/ ) {
            $query_cover      = sprintf "%d", $1;
            $query_cover_perc = sprintf "%d", $2;
        }
        my ( undef, $score ) = split ":", $e[4];
        my ( undef, $rank )  = split ":", $e[5];

        my @pid;
        my @ranges;
        my $hsps = 0;
        my %pos_query;
        my %pos_target;

        while ( ( my $line = <$fh> ) =~ /^HSP/ ) {
            $hsps++;
            if ( $line =~ /pid: (\d+\.?\d+?)$/ ) {
                push @pid, $1;
            }
            if ( $line =~ /query:\((\d+)-(\d+)\)/ ) {
                for ( $1 .. $2 ) {
                    $pos_query{$_}++;
                }
            }
            if ( $line =~ /HSP_ID.+?:\((\d+)-(\d+)\)/ ) {
                for ( $1 .. $2 ) {
                    $pos_target{$_}++;
                }
            }
        }

        my $overlapping = grep { $_ > 1 } values %pos_query;
        my $score_pen = 0;
        if ( $score_strat eq "half" ) {
            $score_pen = sprintf( "%.5f", ( ( $score / $query_cover ) / 2 ) );
        } elsif ( $score_strat eq "full" ) {
            $score_pen = sprintf( "%.5f", ( ( $score / $query_cover ) ) );
        } elsif ( $score_strat eq "extrem" ) {
            $score_pen = 1;
        }
        my $reduce_score_by = $overlapping * $score_pen;
        $score -= $reduce_score_by;
        $score = sprintf( "%.3f", $score );

        my $query_covered  = scalar keys %pos_query;
        my $target_covered = scalar keys %pos_target;

        my $query_length   = $length_query{ $e[0] };
        my $target_length  = $length_target{$target};
        my $query_cov_per  = sprintf "%.2f", $query_covered / $query_length * 100;
        my $target_cov_per = sprintf "%.2f", $target_covered / $target_length * 100;

        my $strand = $e[2] eq "+" ? 1 : -1;

        my $mean_id = sprintf( "%.2f", mean(@pid) );
        print $out join "\t",
          (
            $e[0],                     $target,
            $hsps,                     0,
            $score,                    $score,
            $rank,                     $strand,
            getMin( keys %pos_query ), getMax( keys %pos_query ),
            $target_s,                 $target_e,
            $query_covered,            $target_covered,
            $query_cov_per,            $target_cov_per,
            $mean_id
          );
        print "\n";
    }
}
close $fh;
close $out;

sub mean {
    my $sum = 0;
    for (@_) { $sum += $_; }
    return $sum / @_;
}

1;

__END__


