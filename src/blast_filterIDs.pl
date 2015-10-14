#!/usr/bin/perl 

use strict;
use warnings;

use FindBin;
use lib "$FindBin::Bin/../lib/perl";

use Getopt::Long;
use File::Basename;
use Data::Dumper;

use Bio::SearchIO;
use Bio::SearchIO::Writer::HSPTableWriter;
use Bio::Search::Hit::GenericHit;

use Pod::Usage;

my ($file_blast, $file_ids, $output, $type, $seq);
my $help = 0;
my $man  = 0;
GetOptions(
    "blast=s"  => \$file_blast,
    "ids=s"    => \$file_ids,
    "output=s" => \$output,
    "type=s"   => \$type,
    "seq=s"    => \$seq,
    'help|?'   => \$help,
    man        => \$man
) or pod2usage(1);

pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;

unless (-e $file_blast && -e $file_ids) {
    pod2usage(1);
    die "File not found!\n";
}

my %map;
open FH, "<", $file_ids or die "Can't open $file_ids\n";
while (<FH>) {
    chomp;
    $map{$_} = 1;
}
close FH;

my $index = ($seq eq "query") ? 0 : 1;

open OUT, ">", $output     or die $!;
open FH,  "<", $file_blast or die $!;
while (<FH>) {
    chomp;
    my @e = split "\t";
    my $found = (exists $map{$e[$index]}) ? 1 : 0;

    next if ($type eq "extract" && !$found);
    next if ($type eq "remove"  && $found);

    print OUT $_ . "\n";
}
close FH;
close OUT;

1;

__END__

=head1 NAME

=head1 SYNOPSIS

Extract or remove blast result by subject or query ID.

=head1 OPTIONS

=over 8

=item B<-blast>

BLAST-result (tabular format)

=item B<-ids>

File containing IDs to remove|extract. One ID per line without '>'.

=item B<-output>

Path to resulting BLAST-file.

=item B<-type>

Either 'extract' or 'remove'. Extract keeps only specified IDs.

=item B<-seq>

Search "query" or "subject" for IDs.
    
=back

=cut
