#!/usr/bin/env perl
#
#       AUTHOR: Martin Bens, bensmartin@gmail.com
# ORGANIZATION: FLI Jena
#      CREATED: 03/15/2013 11:16:11
#

use strict;
use warnings;

use Getopt::Long;
use File::Basename;
use Data::Dumper;
use Pod::Usage;

use 5.010;

use constant QUERY_ID   => 0;
use constant TARGET_ID  => 1;
use constant EVALUE     => 3;
use constant IDENTITY   => 16;
use constant QUERY_COV  => 14;
use constant TARGET_COV => 15;

my $file_blast;
my $file_length_query;
my $file_length_subj;

my $thr_evalue       ;
my $thr_identity     ;
my $thr_length_query ;
my $thr_length_subj  ;
my $thr_cov_query    ;
my $thr_cov_subj     ;

my $help = 0;
GetOptions(
    "blast=s"        => \$file_blast,
    "query-index=s"  => \$file_length_query,
    "subj-index=s"   => \$file_length_subj,
    "identity|i=f"           => \$thr_identity,
    "evalue=f"       => \$thr_evalue,
    "query-coverage|qc=f"    => \$thr_cov_query,
    "subj-coverage|sc=f"     => \$thr_cov_subj,
    "query-length|ql=f" => \$thr_length_query,
    "subj-length|sl=f"  => \$thr_length_subj,
    'help|?'         => \$help,
) or pod2usage(-verbose => 2);

pod2usage(-verbose => 2) if $help;

sub getLengthMap {
    my $file = shift;
    my %map;
    open FH, "<", $file or die $!;
    while (<FH>) {
        chomp;
        my @e = split "\t";
        $map{$e[0]} = $e[1];
    }
    close FH;

    return \%map;
}

die "Blast file not found! ($file_blast)\n" unless (-e $file_blast);

# index lenghts for query and subject if required
my ($map_length_query, $map_length_subj);
if (defined $thr_length_query) {
    die "Index file for query id not specified!"
      unless (defined $file_length_query);
    die "Index file not found! ($file_length_query)\n"
      unless (-e $file_length_query);
    $map_length_query = getLengthMap($file_length_query);
}
if (defined $thr_length_subj) {
    die "Index file for subject id not specified!"
      unless (defined $file_length_subj);
    die "Index file not found! ($file_length_subj)\n"
      unless (-e $file_length_subj);
    $map_length_subj = getLengthMap($file_length_subj);
}

open my $fh, "<", $file_blast or die $!;
while (<$fh>) {
    my @e = split "\t";

    if (defined $thr_cov_query) {
        next unless ( $e[QUERY_COV] >= $thr_cov_query );
    }
    if (defined $thr_cov_subj) {
        next unless( $e[TARGET_COV] >= $thr_cov_subj);
    }
    if (defined $thr_length_query) {
        next unless ($map_length_query->{$e[QUERY_ID]} >= $thr_length_query);
    }
    if (defined $thr_length_subj) {
        next unless ($map_length_subj->{$e[TARGET_ID]} >= $thr_length_subj);
    }
    if (defined $thr_evalue) {
        next unless ($e[EVALUE] <= $thr_evalue);
    }
    if (defined $thr_identity) {
        next unless ($e[IDENTITY] >= $thr_identity);
    }

    print $_;

}
close $fh;

1;

__END__

=head1 DESCRIPTION

Filters BLAST-Tabular output (WU BLAST) based on identity, alignment length,
covered fraction of subject or query and evalue. Writes filtered results to
STDOUT and warnings to STDERR.

=head1 EXAMPLE

blast_filter.pl -blast query2subj.csv -id 80.0 -subj-cov 90.0 -subj-index
subj.fai -subj-length 200 1> filtered.csv

=head1 OPTIONS

Only specified filtering options will be used. Comparison is performed with
greater (less) than or equal.

=over 8

=item B<-blast>

WU-Blast output in tabular format.

=item B<-query-index>

Tabular file with id and length in first 2 columns. Used to compute coverage of
query.

=item B<-subj-index>

Tabular file with id and length in first 2 columns. Used to compute coverage of
subject.

=item B<-identity>

Minimum percentage identity.

=item B<-evalue>

Maximum e-value.

=item B<-query-coverage>

Minimum coverage of query sequence (percentage).

=item B<-subj-coverage>

Minimum coverage of subject sequence (percentage).

=item B<-query-length>

Minimum length of query hit.

=item B<-subj-length>

Minimum length of subject hit.

=cut
