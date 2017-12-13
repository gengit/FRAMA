#!/usr/bin/env perl

=pod

=head1 INFO

Martin Bens, bensmartin@gmail.com
2014-11-11

=head1 DESCRIPTION


=head1 OPTIONS

=over 8

=item B<-sequences>

Input data.

=item B<-database>

default: UniVec_Core

=item B<-cpu>

blast parameter

=item B<-outstem>


=back

=cut

use strict;
use warnings;

use Getopt::Long;
use File::Basename;
use Data::Dumper;
use Pod::Usage;
use Bio::SearchIO;

pod2usage( -message => "\n\tNo arguments\n", -verbose => 1 ) if ( @ARGV == 0 );

my $man  = 0;
my $help = 0;

my $seq_file;
my $database = "UniVec_Core";
my $cpus     = 1;
my $ncbi     = 0;
my $outstem;
GetOptions(
    'sequences=s' => \$seq_file,
    'database=s'  => \$database,
    'cpu=i'       => \$cpus,
    'ncbi=i'      => \$ncbi,
    'outstem=s'   => \$outstem,
    'help|?'      => \$help,
    'man|m'       => \$man
) or pod2usage(1);

pod2usage(1) if $help;
pod2usage( -verbose => 2 ) if $man;

my ( $blast_file, $report_file );
unless ($outstem) {
    ( $blast_file  = $seq_file ) =~ s/\.[^\.]+$/.vectorblast/;
    ( $report_file = $seq_file ) =~ s/\.[^\.]+$/.csv/;
} else {
    $blast_file  = "$outstem.vectorblast";
    $report_file = "$outstem.csv";
}

if ($ncbi) {
    unless (-e "$database.nin") {
        my $command = "formatdb -p F -i $database";
        system($command);
    }
    my $command = "blastall -i $seq_file -d $database -a 8 -p blastn -q -5 -G 3 -E 3 -F \"m D\" -e 700 -Y 1.75e12 -o $blast_file";
    system($command);
} else {
    unless (-e "$database.xnd") {
        my $command = "xdformat -n $database";
        system($command);
    }
    my $command = "blastn $database $seq_file -cpus $cpus -M 1 -N -5 -Q 3 -R 3 -wordmask=dust lcmask -E 700 -Y 1.75e12 -o $blast_file &> $blast_file.log";
    system($command);
}

my $blast_stream =
  Bio::SearchIO->new( -format => "blast", -file => "$blast_file" );

open( REP, ">", "$report_file" );
while ( my $result = $blast_stream->next_result() ) {
    my $qname   = $result->query_name;
    my $qlength = $result->query_length;
    my @matches = ();

    while ( my $hit = $result->next_hit() ) {

        while ( my $hsp = $hit->next_hsp() ) {
            my $match_pos =
              (      $hsp->start('query') <= 25
                  || $hsp->end('query') >= $qlength - 25 + 1 )
              ? "terminal"
              : "internal";
            my $match_type = "";
            my $seen       = 0;
            foreach my $m (@matches) {
                if (   $hsp->start('query') >= $m->[0]
                    && $hsp->end('query') <= $m->[1] )
                {
                    $seen = 1;
                    last;
                }
            }
            next if ($seen);

            if ( $match_pos eq "terminal" ) {
                $match_type = (
                    $hsp->score >= 24 ? "strong"
                    : ( $hsp->score >= 19 ? "moderate"
                        : ( $hsp->score >= 16 ? "weak" : "" ) )
                );
            } elsif ( $match_pos eq "internal" ) {
                $match_type = (
                    $hsp->score >= 30 ? "strong"
                    : ( $hsp->score >= 25 ? "moderate"
                        : ( $hsp->score >= 23 ? "weak" : "" ) )
                );
            } else {

            }

            if ($match_type) {
                print REP join( "\t",
                    $qname,               $qlength,
                    $match_pos,           $match_type,
                    $hit->name,           $hit->description,
                    $hsp->start('query'), $hsp->end('query'),
                    $hsp->score() )
                  . "\n";
                push( @matches,
                    [ $hsp->start('query') - 10, $hsp->end('query') + 10 ] );
            }

            last;
        }
    }
}
close(REP);

